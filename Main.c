/************************************************************
			Main.c
************************************************************/

#include"Prop_design.h"

struct loop_counter total_lp,invis_lp;

extern int Set_Design_Elements(struct profile_spec *const spec_PT,struct fixed_points *const r_PT);
extern void Optimize_Gammas(const double RLX_FCTR,double gamma[],const struct fixed_points *const R_PT,struct induced *const id_PT,struct blade_section_data bs_dat[],struct profile_spec *const spec_PT);
extern void Output(const double RLX_FCTR,const struct fixed_points *const R_PT,const struct induced *const ID_PT,const struct blade_section_data BS_DAT[],const double GAM[],const struct profile_spec *const SPEC_PT);
extern void Off_Design_Performance(struct fixed_points *const r_PT,struct blade_section_data bs_dat[],struct induced *const id_PT,double gamma[],struct profile_spec *const spec_PT);
extern void Geometry_design(struct fixed_points *const r_PT,struct profile_spec *const spec_PT);

void Design_procedure(struct fixed_points *const r_PT,struct blade_section_data bs_dat[],struct induced *const id_PT,double gamma[],struct profile_spec *const spec_PT);

struct resultf File;				//結果出力ファイル群
struct airfoil_data Profile_air;	//空力設計部翼型データ
struct airfoil_data Root_air;		//root
char dummy[10];

int main(void)
{
	struct profile_spec spec;
	struct profile_spec *const spec_PT = &spec;
	struct fixed_points r;
	struct fixed_points *const r_PT = &r;	
	struct induced id;
	struct induced *const id_PT = &id;	
	struct blade_section_data bs_dat[JX+1];
	double gamma[JX+1];
	int y_n;
	
	printf("    >=======================================================<\n");
	printf("    >      OPTIMIZATION OF PROPELLERS                       <\n");
	printf("    >              USING HELICOIDAL VORTEX MODEL            <\n");
	printf("    >=======================================================<\n");

	printf("\n計算モード？\n%d:Profile設計から\n%d:Root設計のみ\n%d:性能計算のみ\n>",DESIGN,ROOTDESIGN,OFFDESIGN);
	scanf("%d",&calc_mode);	gets(dummy);

	switch(calc_mode){
	case DESIGN:
		Design_procedure(r_PT,bs_dat,id_PT,gamma,spec_PT);
		Geometry_design(r_PT,spec_PT);
		printf("\n性能計算も行いますか？[y/n]  >");	y_n = getchar();	gets(dummy);
		if(y_n == 'y'){
			calc_mode = OFFDESIGN;
			Off_Design_Performance(r_PT,bs_dat,id_PT,gamma,spec_PT);
		}
		break;
	
	case ROOTDESIGN:
		Read_spec(spec_PT,r_PT);
		Geometry_design(r_PT,spec_PT);
		break;

	case OFFDESIGN:
		Read_spec(spec_PT,r_PT);
		Off_Design_Performance(r_PT,bs_dat,id_PT,gamma,spec_PT);
		break;
	
	default:;
	}

	printf("\n終わり\n\n");
	return 0;
}

void Design_procedure(struct fixed_points *const r_PT,
					  struct blade_section_data bs_dat[],
					  struct induced *const id_PT,
					  double gamma[],
					  struct profile_spec *const spec_PT)
{
	time_t start_time, watch;
	const double Relaxation_factor = 1.2;

	total_lp.relax = 0;		total_lp.lag = 0;
	
	while( Set_Design_Elements(spec_PT,r_PT) == FAILURE );
	
	printf("\nブレード上の循環分布最適化をはじめます\n");
	//初期値設定
	Change_Gammas(INITIALIZE,r_PT,NULL,0,gamma,id_PT,bs_dat,spec_PT);
	
	time(&start_time);	
	//inviscid model
	visc = INVISCID;
		Optimize_Gammas(Relaxation_factor,gamma,r_PT,id_PT,bs_dat,spec_PT);
		time(&watch);	total_lp.time_rec = difftime(watch,start_time);
		invis_lp = total_lp;
	
	//viscous correction
	visc = VISCID;
		Get_blade_section_data(bs_dat,gamma,id_PT,spec_PT);
		Optimize_Gammas(Relaxation_factor,gamma,r_PT,id_PT,bs_dat,spec_PT);
		time(&watch);	total_lp.time_rec = difftime(watch,start_time);

	Evaluate(r_PT,id_PT,gamma,bs_dat,spec_PT);

	//出力
	Output(Relaxation_factor,r_PT,id_PT,bs_dat,gamma,spec_PT);
	printf("Profile設計終了　結果は%sに出力しました\n",File.aero);

}
	
