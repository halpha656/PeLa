/*******************************************************************
		Change_Gammas.c
		初期化も含め設計中gamma[]を変更するのは全部ここでやる
		gamma[]に依存する恒久的変数
		（chord,誘導速度,渦面移動速度,局所空力係数）も計算しなおす
********************************************************************/

#include"Prop_design.h"

void Change_Gammas(enum gam_change_mode gcmode,
				   const struct fixed_points *const R_PT,
				   const double d_gam[JX+1],
				   const double KAPPA,
				   double gamma[],
				   struct induced *const id_PT,
				   struct blade_section_data bs_dat[],
				   struct profile_spec *const spec_PT)
{
	int j;
	
	switch(gcmode){		
		
	case INITIALIZE:	//初期化
	default:
		for(j=1; j<=JX; j++){
			gamma[j] = 0.0;
			spec_PT->chord[j] = 0.5;		//あんま意味無い
			bs_dat[j].c.l_opt = 0.64;
		}
		id_PT->disp_v_half = 0.0;
		break;
		
	case RELAXATION:	//緩和法
		for(j=2; j<=JX-1; j++){
			gamma[j] = gamma[j] + d_gam[j];
		}
		++total_lp.relax;
		break;
		
	case LAGLANGE:		//全体をκ倍
		for(j=2; j<=JX-1; j++)	gamma[j] = gamma[j] * KAPPA;
		++total_lp.lag;
	}

	//gammaによって誘導されるもの諸々を計算
	Biot_Savart(id_PT, R_PT, gamma);
	
	//空力係数etc  非粘性では使わないのでとばす
	if(visc == VISCID)	Get_blade_section_data(bs_dat, gamma, id_PT,spec_PT);

	//chord更新
	for(j=2; j<=JX-1; j++){
		spec_PT->chord[j] = 2 * gamma[j] / ( id_PT->v.q[j] * bs_dat[j].c.l_opt );
		printf("gamma[%d] = %g\tchord[%d] = %g\n",j,gamma[j],j,spec_PT->chord[j]);
	}

}
