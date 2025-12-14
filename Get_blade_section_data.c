/***************************************************************
			Get_blade_section_data.c
			翼素のデータ(structure blade_section_data)
			を求める
****************************************************************/
#include"Prop_design.h"

//線形補間
double Linear_Polation(const double x, const double xref1, const double xref2,
					   const double yref1, const double yref2)
{
	double y;
	y = ( yref1 * ( x - xref2 ) - yref2 * ( x - xref1 ) ) / ( xref1 - xref2 );
	return(y);
}

//局所空力係数（Cdは０～２次のClの係数）、レイノルズ数を計算
void Get_blade_section_data(struct blade_section_data bs_dat[],
							const double GAM[],
							const struct induced *const ID_PT,
							const struct profile_spec *const SPEC_PT)
{	
	int j, m, ii;
	int clx10;
	double gradj, cl0j;
	double clm[3], cdm[3];
	double a,b;

	for(j=2; j<=JX-1; j++){
		
		//翼素のRe
		bs_dat[j].Re_num = SPEC_PT->chord[j] * ID_PT->v.q[j] * Re_prop;
			
		//翼型データの中からReが近い2つ(m番目とm-1番目)を選ぶ	//**なんかあやしい
		for(m=2; m<=Profile_air.num_crvs - 1; m++){
			if(bs_dat[j].Re_num < Profile_air.crvs[m].Re)	break;
			//最後までひっかからなかったら m = num_pola になって出てくる
		}
		
		switch(calc_mode){
		case DESIGN:
			//翼素のCl
			bs_dat[j].c.l = 2 * GAM[j] / ( ID_PT->v.q[j] * SPEC_PT->chord[j] );
		
			//レイノルズ数に応じてCl_optを内挿
			bs_dat[j].c.l_opt = Linear_Polation(bs_dat[j].Re_num,Profile_air.crvs[m].Re, Profile_air.crvs[m-1].Re,
												Profile_air.crvs[m].cl_design, Profile_air.crvs[m-1].cl_design);
			if(bs_dat[j].c.l_opt <= 0.0)	bs_dat[j].c.l_opt = 0.1;//**
			break;

		case OFFDESIGN:
			//迎角α
			Alpha_rad[j] = SPEC_PT->beta[j] - ID_PT->phi[j] + Pitch_rad;
			
			//cl = gradj * α + cl0j
			//データから係数を内挿してαで線形近似	//**要stallの処理
			gradj = Linear_Polation(bs_dat[j].Re_num,Profile_air.crvs[m].Re, Profile_air.crvs[m-1].Re,
									Profile_air.crvs[m].grad_cl_al, Profile_air.crvs[m-1].grad_cl_al);
			cl0j =  Linear_Polation(bs_dat[j].Re_num,Profile_air.crvs[m].Re, Profile_air.crvs[m-1].Re,
									Profile_air.crvs[m].cl_al0, Profile_air.crvs[m-1].cl_al0);
			bs_dat[j].c.l = gradj * Alpha_rad[j] + cl0j;
			break;

		default:
			printf("Error!\nGet_blade_section_data\ncalc_mode = %d\n",calc_mode);
			exit(EXIT_FAILURE);
		}	
	
		//cl*10を四捨五入
		clx10 = (int)(bs_dat[j].c.l * 10.0 + 0.5);
		
		//clがデータ範囲外にはみ出たらはじっこに強制--->Cdは外挿になる
		if( clx10 >= Profile_air.polar_size ) clx10 = Profile_air.polar_size - 1;		
		if( clx10 <= 0 )	clx10 = 1;

		for(ii=0; ii<=2; ii++){
			clm[ii] = (double)(clx10 - 1 + ii) / 10.0;		
			cdm[ii] = Linear_Polation(bs_dat[j].Re_num, Profile_air.crvs[m].Re, Profile_air.crvs[m-1].Re,
								Profile_air.crvs[m].cd[clx10-1+ii], Profile_air.crvs[m-1].cd[clx10-1+ii]);
		}
	
		//CdをClで二次近似	係数を求める
		//cd = cdm[0] + cdm[1] * cl + cdm[2] * cl^2
		a = ( cdm[1] - cdm[0] ) / ( clm[1] - clm[0] );
		b = ( (cdm[2] - cdm[1])*(clm[1] - clm[0]) - (cdm[1] - cdm[0])*(clm[2] - clm[1]) )
			/ ( (clm[2] - clm[1])*(clm[2] - clm[0])*(clm[1] - clm[0]) );
		
		bs_dat[j].c.dm[0] = cdm[0] - a * clm[0] + b * clm[0] * clm[1];
		bs_dat[j].c.dm[1] = a - b * (clm[0] + clm[1]);
		bs_dat[j].c.dm[2] = b;

		bs_dat[j].c.d = bs_dat[j].c.dm[0] 
						+ bs_dat[j].c.dm[1] * bs_dat[j].c.l 
						+ bs_dat[j].c.dm[2] * bs_dat[j].c.l * bs_dat[j].c.l;
		
	}

}
