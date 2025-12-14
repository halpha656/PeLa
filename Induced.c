/********************************************************
		Induced.c
		循環分布から誘導速度etc（structure induced）
		に値を与える
********************************************************/

#include"Prop_design.h"

//**一番時間がかかる箇所	
void Factor(const struct fixed_points *const R_PT,struct induced *const id_PT)
{
	int    i, j, k;
	double y_v[IX+1], z_v[IX+1];
	double y_v_ave, z_v_ave;
	double r_v0, r_v1;
	double adv_var;

	//a[j][k],b[j][k]を求める（ 1 はもう一枚のブレードからの影響）
	for(j=1; j<=JX-1; j++){
			
		for(k=1; k<=JX-1; k++){
				
			id_PT->fc.a0[j][k] = 0.0;	id_PT->fc.a1[j][k] = 0.0;
			id_PT->fc.b0[j][k] = 0.0;	id_PT->fc.b1[j][k] = 0.0;
			
			for(i=1; i<=IX; i++){
				adv_var = Adv * ( 1.0 + id_PT->disp_v_half * R_PT->accel_mdl[i]);
				y_v[i] = R_PT->y[j] * cos( R_PT->x_v[i] / adv_var );
				z_v[i] = R_PT->y[j] * sin( R_PT->x_v[i] / adv_var );
			}
			for(i=2; i<=IX; i++){			
				
				y_v_ave = (y_v[i] + y_v[i-1])/2.0;
				z_v_ave = (z_v[i] + z_v[i-1])/2.0;

				r_v0 = sqrt(
							(R_PT->x_v[i-1] + R_PT->x_v[i]) * (R_PT->x_v[i-1] + R_PT->x_v[i])/4.0 
							+ ( R_PT->y_p[k] - y_v_ave ) * ( R_PT->y_p[k] - y_v_ave ) 
							+ z_v_ave * z_v_ave 
							);
				r_v1 = sqrt(
							(R_PT->x_v[i-1] + R_PT->x_v[i]) * (R_PT->x_v[i-1] + R_PT->x_v[i])/4.0 
							+ ( R_PT->y_p[k] + y_v_ave ) * ( R_PT->y_p[k] + y_v_ave ) 
							+ z_v_ave * z_v_ave
							);
				
				id_PT->fc.a0[j][k] = id_PT->fc.a0[j][k] - ( ( R_PT->y_p[k] - y_v_ave) * (R_PT->x_v[i] - R_PT->x_v[i-1]) + (R_PT->x_v[i-1] + R_PT->x_v[i])/2.0 * (y_v[i] - y_v[i-1]) ) / ( r_v0 * r_v0 * r_v0 );
				id_PT->fc.a1[j][k] = id_PT->fc.a1[j][k] - ( ( R_PT->y_p[k] + y_v_ave) * (R_PT->x_v[i] - R_PT->x_v[i-1]) + (R_PT->x_v[i-1] + R_PT->x_v[i])/2.0 * (y_v[i-1] - y_v[i]) ) / ( r_v1 * r_v1 * r_v1 );
				id_PT->fc.b0[j][k] = id_PT->fc.b0[j][k] + ( ( R_PT->y_p[k] - y_v_ave) * (z_v[i] - z_v[i-1]) + z_v_ave * (y_v[i] - y_v[i-1]) ) / ( r_v0 * r_v0 * r_v0 );
				id_PT->fc.b1[j][k] = id_PT->fc.b1[j][k] + ( ( R_PT->y_p[k] + y_v_ave) * (-z_v[i] + z_v[i-1]) - z_v_ave * (y_v[i-1] - y_v[i]) ) / ( r_v1 * r_v1 * r_v1 );
			}		

			id_PT->fc.a0[j][k] = id_PT->fc.a0[j][k] / ( 4.0 * PI );		id_PT->fc.b0[j][k] = id_PT->fc.b0[j][k] / ( 4.0 * PI );
			id_PT->fc.a1[j][k] = id_PT->fc.a1[j][k] / ( 4.0 * PI );		id_PT->fc.b1[j][k] = id_PT->fc.b1[j][k] / ( 4.0 * PI );
		}
	}
}

//誘導速度w,u、局所流入速度q、流入角φを変更
void Biot_Savart(struct induced *const id_PT, const struct fixed_points *const R_PT, const double GAM[])
{
	int j, k;

	Factor(R_PT,id_PT);
	
	for(k=1; k<=JX-1; k++){
		
		id_PT->v.w[k] = 0.0;		id_PT->v.u[k] = 0.0;
		
		for(j=1; j<=JX-1; j++){
			id_PT->v.w[k] = id_PT->v.w[k] + ( GAM[j+1] - GAM[j] ) * (id_PT->fc.a0[j][k] + id_PT->fc.a1[j][k]) ;
			id_PT->v.u[k] = id_PT->v.u[k] + ( GAM[j+1] - GAM[j] ) * (id_PT->fc.b0[j][k] + id_PT->fc.b1[j][k]) ;
		}
		//合成流入速度[-]		
		id_PT->v.q[k] = sqrt( 
							(1 + id_PT->v.u[k]) * (1 + id_PT->v.u[k]) 
							+ ( R_PT->y_p[k] / Adv + id_PT->v.w[k] ) * ( R_PT->y_p[k] / Adv + id_PT->v.w[k] )
							);
		//流入角[rad]
		id_PT->phi[k] = atan( (1+id_PT->v.u[k])/(R_PT->y_p[k]/Adv + id_PT->v.w[k]) );
	}
	
	//**displacement velocity　代表として７０％あたりを使う  適当
	//**始めからやると収束しないので非粘性の時は無視(ゼロ)
	if(visc == VISCID)	id_PT->disp_v_half = id_PT->v.u[J70] - id_PT->v.w[J70] * tan(id_PT->phi[J70]);

}

/**********************************************************
		性能計算用	誘導速度計算
**********************************************************/

void Momentum_loss_function(const struct fixed_points *const R_PT, const double phi_tip,
							struct induced *const id_PT)
{
	int j;
	double f;

	for(j=1; j<=JX-1; j++){
		f = ( 1.0 - R_PT->y_p[j] ) / sin(phi_tip);
		id_PT->loss_func[j] = atan( sqrt( exp(2.0 * f) - 1 ) ) * 2.0 / PI;
	}

}

void Adkins_Liebeck(struct induced *const id_PT, 
					const double LOC_SOLID[],
					const struct blade_section_data BS_DAT[],
					const struct fixed_points *const R_PT)
{
	double coef_x;
	double coef_z;
	double kx;
	double kz;
	double bs_cd[JX+1];		//blade section drag coefficient
	int j;

	//速度減少関数を求める
	Momentum_loss_function(R_PT,id_PT->phi[JX-1],id_PT);

	//誘導速度を求める
	for(j=2; j<=JX-1; j++){
		bs_cd[j] = BS_DAT[j].c.dm[0] 
					+ BS_DAT[j].c.dm[1] * BS_DAT[j].c.l 
					+ BS_DAT[j].c.dm[2] * BS_DAT[j].c.l * BS_DAT[j].c.l;

		coef_x = - BS_DAT[j].c.l * cos(id_PT->phi[j]) + bs_cd[j] * sin(id_PT->phi[j]);
		coef_z = BS_DAT[j].c.l * sin(id_PT->phi[j]) + bs_cd[j] * cos(id_PT->phi[j]);

		kx = - coef_x / ( 4.0 * sin(id_PT->phi[j]) * sin(id_PT->phi[j]) );
		kz =   coef_z / ( 4.0 * sin(id_PT->phi[j]) * cos(id_PT->phi[j]) );

		id_PT->v.u[j] =   LOC_SOLID[j] * kx / ( id_PT->loss_func[j] - LOC_SOLID[j] * kx );
		id_PT->v.w[j] = - LOC_SOLID[j] * kz / ( id_PT->loss_func[j] + LOC_SOLID[j] * kz );

		//合成流入速度[-]		
		id_PT->v.q[j] = sqrt( 
							(1 + id_PT->v.u[j]) * (1 + id_PT->v.u[j]) 
							+ ( R_PT->y_p[j] / Adv + id_PT->v.w[j] ) * ( R_PT->y_p[j] / Adv + id_PT->v.w[j] )
							);
		//流入角[rad]
		id_PT->phi[j] = atan( (1+id_PT->v.u[j])/(R_PT->y_p[j]/Adv + id_PT->v.w[j]) );
	}
}

		