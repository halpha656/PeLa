/**********************************************************
	Off_Design_Performance.c
	設計点外の性能計算
**********************************************************/

#include"Prop_design.h"

extern void Set_Off_Design_Condition(const struct profile_spec *const SPEC_PT);
extern void Adkins_Liebeck(struct induced *const id_PT,const double LOC_SOLID[],const struct blade_section_data BS_DAT[],const struct fixed_points *const R_PT);
extern void Output2(const struct fixed_points *const R_PT,const struct induced *const ID_PT,const struct blade_section_data BS_DAT[],const double GAM[]);

void Read_spec(struct profile_spec *const spec_PT,struct fixed_points *const r_PT);
void State_of_Operation(const struct fixed_points *const R_PT,struct blade_section_data bs_dat[],struct induced *const id_PT,double gamma[],const struct profile_spec *const SPEC_PT);
void Local_solidity(const struct fixed_points *const R_PT,const struct profile_spec *const SPEC_PT,double loc_solid[]);

void Off_Design_Performance(struct fixed_points *const r_PT,
							struct blade_section_data bs_dat[],
							struct induced *const id_PT,
							double gamma[],
							struct profile_spec *const spec_PT)
{
	int y_n;

	visc = VISCID;
	do{
		Set_Off_Design_Condition(spec_PT);
		State_of_Operation(r_PT,bs_dat,id_PT,gamma,spec_PT);
		Evaluate(r_PT,id_PT,gamma,bs_dat,spec_PT);
		Output2(r_PT,id_PT,bs_dat,gamma);
		printf("\n性能計算を続けますか？[y/n]  >");	y_n = getchar();	gets(dummy);
	}while(y_n == 'y');

}

void State_of_Operation(const struct fixed_points *const R_PT,
						struct blade_section_data bs_dat[],
						struct induced *const id_PT,
						double gamma[],
						const struct profile_spec *const SPEC_PT)
{
	int j;
	int n = 0;
	double loc_solid[JX+1];

	Local_solidity(R_PT,SPEC_PT,loc_solid);
	
	//誘導速度初期値設定
	for(j=1; j<=JX-1; j++){
		id_PT->v.u[j] = 0.0;	id_PT->v.w[j] = 0.0;
		id_PT->v.q[j] = sqrt( 1 + R_PT->y_p[j] * R_PT->y_p[j] / (Adv * Adv ));
		id_PT->phi[j] = atan( Adv / R_PT->y_p[j] );
	}
	id_PT->disp_v_half = 0.0;

	while(n < 10){	//適当
		++n;

		Get_blade_section_data(bs_dat,NULL,id_PT,SPEC_PT);
		
		//gamma更新
		for(j=2; j<=JX-1; j++)	gamma[j] = id_PT->v.q[j] * SPEC_PT->chord[j] * bs_dat[j].c.l / 2.0;

		//誘導速度を計算
		Adkins_Liebeck(id_PT,loc_solid,bs_dat,R_PT);
	}
}

void Local_solidity(const struct fixed_points *const R_PT,
					const struct profile_spec *const SPEC_PT,
					double loc_solid[])
{
	int j;
	for(j=1; j<=JX; j++)	loc_solid[j] = SPEC_PT->chord[j] / ( R_PT->y_p[j] * PI );	
}