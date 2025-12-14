/************************************************************
		Optimeze_Gamma_Distribution.c
		緩和法を用いた循環分布最適化
************************************************************/
#include"Prop_design.h"

//翼素の循環で偏微分された係数(Relaxtion form)
static struct pd_Coefs{
	double dg;		//drag(negative thrust)coefficient　の偏微分
	double tq;		//torque coefficientの偏微分
	double dg_gam;
	double tq_gam;
};
void Relaxation(const double RLX_FCTR,double gamma[],const struct fixed_points *const R_PT,struct induced *const id_PT,struct blade_section_data bs_dat[],const double LAMDA,struct profile_spec *const spec_PT);
struct pd_Coefs Partial_Differential(const int J,const struct blade_section_data BS_DATJ,const double GAM[],const struct fixed_points *const R_PT,const struct induced *const ID_PT,const struct profile_spec *const SPEC_PT);

//循環分布を最適化
void Optimize_Gammas(const double RLX_FCTR,
					 double gamma[],
					 const struct fixed_points *const R_PT,
					 struct induced *const id_PT,
					 struct blade_section_data bs_dat[],
					 struct profile_spec *const spec_PT)
{
	double kappa = 0.0;			//Γを変える割合
	static double lamda = 0.0;	//トルクにかけるウエイト
	struct Coefs co;

	while(fabs(1.0 - kappa) > CONVERGE_LAG){

		Relaxation(RLX_FCTR,gamma,R_PT,id_PT,bs_dat,lamda,spec_PT);
		
		co = Coefficients(bs_dat,gamma,R_PT,id_PT,spec_PT);
		
		kappa = ( sqrt( co.tq[1] * co.tq[1] + 4 * co.tq[2] * ( CT_target - co.tq[0] ) ) - co.tq[1] ) / ( 2 * co.tq[2] );
		Change_Gammas(LAGLANGE,R_PT,NULL,kappa,gamma,id_PT,bs_dat,spec_PT);
		printf("kappa = %g\nVISC = %d\n",kappa,visc);

		lamda = (-1.0 ) * ( co.dg[1] + 2 * kappa * co.dg[2] ) / ( co.tq[1] + 2 * kappa * co.tq[2] );
	}	

}

//緩和法 でΓ更新
void Relaxation(const double RLX_FCTR,
				double gamma[],
				const struct fixed_points *const R_PT,
				struct induced *const id_PT,
				struct blade_section_data bs_dat[],
				const double LAMDA,
				struct profile_spec *const spec_PT)
{
	int j;
	double dgam[JX+1];
	double sum_gamma = 1.0;
	double sum_dgam=1000, sum_dgam0=100000;
	struct pd_Coefs pdC;

	while((sum_dgam/sum_gamma > CONVERGE_RLX) && (sum_dgam < sum_dgam0)){	//収束判定
		
		sum_dgam0 = sum_dgam;	sum_dgam = 0.0;
		sum_gamma = 0.0;
		for(j=2; j<=JX-1; j++){
			pdC = Partial_Differential(j,bs_dat[j],gamma,R_PT,id_PT,spec_PT);

			dgam[j] = - RLX_FCTR * ( pdC.dg + LAMDA * pdC.tq ) / ( pdC.dg_gam + LAMDA * pdC.tq_gam );
			sum_dgam = sum_dgam + fabs(dgam[j]) * R_PT->dy[j];
			sum_gamma = sum_gamma + fabs(gamma[j]) * R_PT->dy[j];
		}

		Change_Gammas(RELAXATION,R_PT,dgam,0,gamma,id_PT,bs_dat,spec_PT);
		printf("sum_gamma = %g\tsum_dgam = %g\tsum_dgam0 = %g\nLAMDA = %g\nvisc = %d\ndisp_v_half = %g\n",
			sum_gamma,sum_dgam,sum_dgam0,LAMDA,visc,id_PT->disp_v_half);
	}
}

//プロペラ全体の抵抗係数、トルク係数を求める（積分）
struct Coefs Coefficients(const struct blade_section_data BS_DAT[],
						  const double GAM[],
						  const struct fixed_points *const R_PT,
						  const struct induced *const ID_PT,
						  const struct profile_spec *const SPEC_PT)
{
	struct Coefs co;
	int k,ii;

	for(ii=0; ii<=2; ii++){
		co.dg[ii] = 0.0;
		co.tq[ii] = 0.0;
	}

	for(k=2; k<=JX-1; k++){
		co.dg[1] = co.dg[1] + GAM[k] * R_PT->y_p[k] / Adv * R_PT->dy[k];
		co.dg[2] = co.dg[2] + GAM[k] * ID_PT->v.w[k] * R_PT->dy[k];
		
		co.tq[2] = co.tq[2] + GAM[k] * ID_PT->v.u[k] * ( R_PT->y_p[k] ) * R_PT->dy[k];
	}

	co.dg[1] = co.dg[1] * (-4.0);
	co.dg[2] = co.dg[2] * (-4.0);
	co.tq[1] = co.dg[1] * (-Adv);
	co.tq[2] = co.tq[2] * 4.0;

	//粘性アリ
	if(visc == VISCID){
		for(ii=0; ii<=2; ii++){
			for(k=2; k<=JX-1; k++){
				co.dg[ii] = co.dg[ii] + 2 * ID_PT->v.q[k] * ( 1 + ID_PT->v.u[k] )
							* BS_DAT[k].c.dm[ii] * pow( BS_DAT[k].c.l, ii ) * SPEC_PT->chord[k] * R_PT->dy[k];
				
				co.tq[ii] = co.tq[ii] + 2 * ID_PT->v.q[k] * ( R_PT->y_p[k]/Adv + ID_PT->v.w[k] )
							* BS_DAT[k].c.dm[ii] * pow( BS_DAT[k].c.l, ii ) * R_PT->y_p[k] * SPEC_PT->chord[k] * R_PT->dy[k];
			}
		}
	}

	co.dgtotal = co.dg[0] + co.dg[1] + co.dg[2];
	co.tqtotal = co.tq[0] + co.tq[1] + co.tq[2];

	return(co);
}

//抵抗係数、トルク係数をJ番目の翼素のΓで変微分(Relaxtion form)
struct pd_Coefs Partial_Differential(const int J,
									 const struct blade_section_data BS_DATJ,
									 const double GAM[],
									 const struct fixed_points *const R_PT,
									 const struct induced *const ID_PT,
									 const struct profile_spec *const SPEC_PT)
{
	struct pd_Coefs pdC;
	
	int k;
	double sum_dg = 0.0;
	double sum_tq = 0.0;

	for(k=2; k<=JX-1; k++){
		sum_dg = sum_dg 
				 + GAM[k] 
					* ( ID_PT->fc.a0[J-1][k] + ID_PT->fc.a1[J-1][k] - ID_PT->fc.a0[J][k] - ID_PT->fc.a1[J][k])
					* R_PT->dy[k];
		
		sum_tq = sum_tq
				 + GAM[k]
					* ( ID_PT->fc.b0[J-1][k] + ID_PT->fc.b1[J-1][k] - ID_PT->fc.b0[J][k] - ID_PT->fc.b1[J][k])
					* R_PT->y_p[k] * R_PT->dy[k]; 
	}

	pdC.dg = - 4.0 * ( R_PT->y_p[J] / Adv + ID_PT->v.w[J] ) * R_PT->dy[J] - 4 * sum_dg ;
	pdC.tq =   4.0 * ( 1 + ID_PT->v.u[J] ) * R_PT->y_p[J] * R_PT->dy[J] + 4 * sum_tq ;

	pdC.dg_gam = -8 * ( ID_PT->fc.a0[J-1][J] - ID_PT->fc.a0[J][J] ) * R_PT->dy[J];
	pdC.tq_gam =  8 * ( ID_PT->fc.b0[J-1][J] - ID_PT->fc.b0[J][J] ) * R_PT->y_p[J] * R_PT->dy[J];

	if(visc == VISCID){
		pdC.dg = pdC.dg 
				+ 4 * ( 1 + ID_PT->v.u[J] ) 
					* ( BS_DATJ.c.dm[1] + 4 * BS_DATJ.c.dm[2] * GAM[J] / ( ID_PT->v.q[J] * SPEC_PT->chord[J] ) ) 
					* R_PT->dy[J];
		
		pdC.tq = pdC.tq 
				+ 4 * ( R_PT->y_p[J] / Adv + ID_PT->v.w[J] ) 
					* ( BS_DATJ.c.dm[1] + 4 * BS_DATJ.c.dm[2] * GAM[J] / ( ID_PT->v.q[J] * SPEC_PT->chord[J] ) )
					* R_PT->y_p[J] * R_PT->dy[J];
	
		pdC.dg_gam = pdC.dg_gam 
					+ 16 * ( 1 + ID_PT->v.u[J] ) * BS_DATJ.c.dm[2] 
						* R_PT->dy[J] / ( ID_PT->v.q[J] * SPEC_PT->chord[J] );
		
		pdC.tq_gam = pdC.tq_gam 
					+ 16 * ( R_PT->y_p[J] / Adv + ID_PT->v.w[J] ) * BS_DATJ.c.dm[2] * R_PT->y_p[J] 
						* R_PT->dy[J] / ( ID_PT->v.q[J] * SPEC_PT->chord[J] );
	}

	return(pdC);
}
