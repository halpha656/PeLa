/***************************************************
		Interface.c
		入出力  有次元
***************************************************/

#include"Prop_design.h"

double Re_prop;				//基準レイノルズ数
double Adv;					//進行率[-]  V/(RΩ)
double CT_target;			//トルク係数[-]

double Alpha_rad[JX+1];		//有効迎角[rad]
double Pitch_rad;			//ピッチ角[rad]

static struct range attack_angle;
static struct range stall_section;
static struct range ngtv_lift_section;

static double RHO = 1.184;				//空気密度[kg/m^3]
static double NU = 1.5545e-5;			//動粘性係数[m^2/s]

static double V_cruise;				//機体巡航速度[m/s]
static double P_input;				//入力パワ[W]
static double Rpm;					//回転数[rpm]
static double Design_al_deg;		//設計有効迎角[deg]
static double Omega;				//回転角速度[rad/s]
static double Pitch_deg;			//ピッチ角[deg]

static double Thrust;				//推力[N]　
static double Eta;					//プロペラ効率[-]

static struct Coefs Co_product;
extern void Write_spec(FILE *const aerof_pt, const struct profile_spec *const SPEC_PT);

int Set_Design_Elements(struct profile_spec *const spec_PT,struct fixed_points *const r_PT)
{
	int y_n;
	double y0;
	char resourcef[100];

	printf("\n設計モード？\n%d:新規設計\n%d:Canopus\n%d:test\n>",NEW,CANOPUS,TEST);
	scanf("%d",&design_mode);	gets(dummy);

	switch(design_mode){
	case NEW:
		printf("命名                  >");	scanf("%s",File.propname);
		printf("巡航速度[m/s]?        >");	scanf("%lf",&V_cruise);
		printf("入力パワ[W]?          >");	scanf("%lf",&P_input);
		printf("回転数[rpm]?          >");	scanf("%lf",&Rpm);
		printf("回転半径[m]?          >");	scanf("%lf",&(spec_PT->radius));
		printf("設計内側半径[m]?      >");	scanf("%lf",&(spec_PT->r0));
		printf("翼型データファイル名？>");	scanf("%s",resourcef);
		printf("有効迎角[deg]?        >");	scanf("%lf",&Design_al_deg);
		gets(dummy);
		
		printf("入力ミス？(ある:y　ない:n)	>");	y_n = getchar();	gets(dummy);
		if(y_n == 'y')	return FAILURE;
		break;

	case CANOPUS:
		strcpy(File.propname,"canopus");
		V_cruise = 7.0;	
		Rpm = 140.0;
		spec_PT->radius = 1.65;
		P_input = 270.0;
		spec_PT->r0 = 0.30;
		strcpy(resourcef,"DAE51ex.txt");
		Design_al_deg = 0.0;
		break;

	case TEST:
	default:
		strcpy(File.propname,"test");
		V_cruise = 7.0;
		P_input = 280;
		Rpm = 140;
		spec_PT->radius = 1.5;
		spec_PT->r0 = 0.2;
		strcpy(resourcef,"FORTEST.txt");
		Design_al_deg = 1.0;
		break;
	
	}

	printf("巡航速度     %.1f[m/s]\n",V_cruise);
	printf("回転数       %.0f[rpm]\n",Rpm);
	printf("回転半径     %.2f[m]\n",spec_PT->radius);
	printf("入力パワ     %.0f[W]\n",P_input);
	printf("設計内側半径 %.2f[m]\n",spec_PT->r0);
	printf("翼型データ %s\n",resourcef);
	printf("有効迎角     %.1f[deg]\n",Design_al_deg);

	y0 = spec_PT->r0 / spec_PT->radius;
	spec_PT->design_al_rad = Design_al_deg * PI / 180.0;

	//設計用翼型データの読み込み
	Profile_air = Read_Airfoil_Data(PROFILE, resourcef,spec_PT);
	if( strcmp(Profile_air.resource,"unknown") == 0 )	return FAILURE;

	//計算点座標を設定
	Set_point(r_PT,spec_PT);
	
	//定数設定
	Re_prop = V_cruise * spec_PT->radius / NU;
	Omega = 2.0 * PI * Rpm / 60.0;
	Adv = V_cruise / ( spec_PT->radius * Omega );
	CT_target = 2.0 * P_input / ( RHO * Omega * pow(spec_PT->radius,3) * V_cruise * V_cruise);

	//ファイル作成
	if(Make_files() != FAILURE)	return 0;
	else return(FAILURE);
	
}

int Make_files(void)
{
	FILE *fp;

	sprintf(File.aero,"%s_aero.csv",File.propname);
	fp = fopen(File.aero,"a");
	if(fp == NULL){
		printf("%sを開けません\n");
		return FAILURE;
	}
	else{
		fclose(fp);
		sprintf(File.geometry,"%s_geom.csv",File.propname);
		sprintf(File.cad,"%s_cad.scr",File.propname);
		return 0;
	}
}

//推力、パワ、効率を計算　次元つける
void Evaluate(const struct fixed_points *const R_PT,
			  const struct induced *const ID_PT,
			  const double GAM[],
			  const struct blade_section_data BS_DAT[],
			  struct profile_spec *const spec_PT)
{		

	int j;

	attack_angle.lower = 10.0;	attack_angle.upper = -10.0;
	stall_section.lower = 10.0;	stall_section.upper = -10.0;
	ngtv_lift_section.lower = 10.0;	ngtv_lift_section.upper = -10.0;

	Co_product = Coefficients(BS_DAT,GAM,R_PT,ID_PT,spec_PT);

	//推力
	Thrust = - 1 / 2.0 * Co_product.dgtotal * RHO * V_cruise * V_cruise 
					* spec_PT->radius * spec_PT->radius;
	//プロペラ効率(有効仕事/入力パワー)
	Eta = - Adv * ( Co_product.dgtotal / Co_product.tqtotal );
	
	switch(calc_mode){
	case DESIGN:
		//取り付け角[rad]
		for(j=2; j<=JX-1; j++)	spec_PT->beta[j] = ID_PT->phi[j] + spec_PT->design_al_rad;
		break;

	case OFFDESIGN:
		//必要な入力パワー
		P_input = 1 / 2.0 * Co_product.tqtotal * RHO * V_cruise * V_cruise 
						* spec_PT->radius * spec_PT->radius * spec_PT->radius * Omega ;
		
		//大体の作動状況　迎角分布・失速域・負揚力域（めやす）
		for(j=2; j<=JX-3; j++){
			if(Alpha_rad[j] >= attack_angle.upper)	attack_angle.upper = Alpha_rad[j];
			if(Alpha_rad[j] <= attack_angle.lower)	attack_angle.lower = Alpha_rad[j];
			if(Alpha_rad[j] >= BS_DAT[j].stall){
				if(R_PT->y_p[j] < stall_section.lower)	stall_section.lower = R_PT->y_p[j];
				if(R_PT->y_p[j] > stall_section.upper)	stall_section.upper = R_PT->y_p[j];
			}
			if(BS_DAT[j].c.l < 0.0){
				if(R_PT->y_p[j] < ngtv_lift_section.lower)	ngtv_lift_section.lower = R_PT->y_p[j];
				if(R_PT->y_p[j] > ngtv_lift_section.upper)	ngtv_lift_section.upper = R_PT->y_p[j];
			}
		}

	default:;
	}


	printf("\n<結果>\n");
	printf("推力\t%g[N]\t%g[kgf]\n",Thrust,Thrust/9.8);
	printf("入力パワ\t%g[W]\n",P_input);
	printf("有効仕事\t%g[W]\n",Thrust * V_cruise);
	printf("プロペラ効率\t%g[-]\n",Eta);
	
	if(calc_mode == OFFDESIGN){
		printf("\n有効迎角\t%g～%g[deg]\n",attack_angle.lower*180/PI,attack_angle.upper*180/PI);
		
		printf("失速領域\t");
		if(stall_section.lower == 10.0)	printf("なし\n");
		else printf("%g～%g[%%]\n",stall_section.lower*100, stall_section.upper*100);
		
		printf("負揚力域\t");
		if(ngtv_lift_section.lower == 10.0)	printf("なし\n");
		else printf("%g～%g[%%]\n",ngtv_lift_section.lower*100, ngtv_lift_section.upper*100);
	}

}

void Set_Off_Design_Condition(const struct profile_spec *const SPEC_PT)
{
	printf("\nフライト条件を設定してください\n");
	printf("機体速度[m/s]? >");	scanf("%lf",&V_cruise);
	printf("回転数[rpm]?   >");	scanf("%lf",&Rpm);
	printf("ピッチ[deg]?   >");	scanf("%lf",&Pitch_deg);
	gets(dummy);

	Pitch_rad = Pitch_deg * PI / 180.0;
	Re_prop = V_cruise * SPEC_PT->radius / NU;
	Omega = 2.0 * PI * Rpm / 60.0;
	Adv = V_cruise / ( SPEC_PT->radius * Omega );
}

void Output(const double RLX_FCTR,
			const struct fixed_points *const R_PT,
			const struct induced *const ID_PT,
			const struct blade_section_data BS_DAT[],
			const double GAM[],
			const struct profile_spec *const SPEC_PT)
{
	FILE *aerof_pt;
	int j;
	
	aerof_pt = fopen(File.aero,"w");
	//計算結果
	fprintf(aerof_pt,"\nρ,%g,[kg/m^3],ν,%g,[m^2/s]\n",RHO,NU);
	fprintf(aerof_pt,"<入力パラメータ>\n");
	fprintf(aerof_pt,"巡航速度,%g,[m/s]\n",V_cruise);
	fprintf(aerof_pt,"回転半径,%g,[m],設計内側半径,%g,[m]\n回転数,%g,[rpm]\n翼型,%s\n有効迎角,%g,[deg]\n進行率,%g,[-]\n",SPEC_PT->radius,SPEC_PT->r0,Rpm,Profile_air.name,Design_al_deg,Adv);
	fprintf(aerof_pt,"入力パワ,%g,[W],Cτ target,%g,[-]\n",P_input,CT_target);
	
	fprintf(aerof_pt,"\n<計算環境>\n");
	fprintf(aerof_pt,"Relaxation_factor,%g,[-]\n",RLX_FCTR);
	fprintf(aerof_pt,"分割,JX,%d,IX,%d,SX,%g,DX_0,%g\n",JX,IX,SX,DX_0);
	fprintf(aerof_pt,"\n収束判定,CONVERGE_RLX,%g,CONVERGE_LAG,%g\n",CONVERGE_RLX,CONVERGE_LAG);

	fprintf(aerof_pt,"\n<結果>\n");
	fprintf(aerof_pt,"推力,%g,[N],%g,[kgf]\n",Thrust,Thrust/9.8);
	fprintf(aerof_pt,"有効仕事,%g,[W]\n",Thrust * V_cruise);
	fprintf(aerof_pt,"プロペラ効率,%g,[-]\n",Eta);
	fprintf(aerof_pt,"後流移動速度 /V [-],%g\n",ID_PT->disp_v_half);
	
	fprintf(aerof_pt,"\n翼素データ\n");
	fprintf(aerof_pt,"y_p/R[-],Γ/(VR)[-],u/V[-],w/V[-],q/V[-],Re[-],chord/R[-],phi[rad],beta[rad],phi[deg],beta[deg],Cl_opt[-],Cl[-],Cd[-],u - w * tanφ\n");
	for(j=2; j<=JX-1; j++){
		fprintf(aerof_pt,"%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n"
			,R_PT->y_p[j],GAM[j],ID_PT->v.u[j],ID_PT->v.w[j],ID_PT->v.q[j],
			BS_DAT[j].Re_num,SPEC_PT->chord[j],ID_PT->phi[j],SPEC_PT->beta[j],ID_PT->phi[j]*180/PI,SPEC_PT->beta[j]*180/PI,
			BS_DAT[j].c.l_opt,BS_DAT[j].c.l,
			BS_DAT[j].c.dm[0] + BS_DAT[j].c.dm[1]*BS_DAT[j].c.l + BS_DAT[j].c.dm[2]*BS_DAT[j].c.l*BS_DAT[j].c.l,
			ID_PT->v.u[j] - ID_PT->v.w[j]*tan(ID_PT->phi[j]));
	}
	
	fprintf(aerof_pt,"\nΓ変更回数,,,,,,時間[s]\n非粘性,Relaxation,%d,Laglange,%d,,%f\n合計,Relaxation,%d,Laglange,%d,,%f\n\n",invis_lp.relax,invis_lp.lag,invis_lp.time_rec,total_lp.relax,total_lp.lag,total_lp.time_rec);
	
	Write_spec(aerof_pt, SPEC_PT);
	
	fclose(aerof_pt);
}

void Output2(const struct fixed_points *const R_PT,
			 const struct induced *const ID_PT,
			 const struct blade_section_data BS_DAT[],
			 const double GAM[])
{
	FILE *aerof_pt;
	int j;
	int y_n;

	aerof_pt = fopen(File.aero,"a");
	fprintf(aerof_pt,"\n\n>>>>>>>>>>>>>Adkins&Liebeckの理論による設計点外性能計算>>>>>>>>>>>>>>>>\n");
	fprintf(aerof_pt,"<フライト条件>\n");
	fprintf(aerof_pt,"機体速度,%g,[m/s]\n回転数,%g,[rpm]\nピッチ角,%g,[deg]\n",V_cruise,Rpm,Pitch_deg);
	
	fprintf(aerof_pt,"<結果>\n");
	fprintf(aerof_pt,"推力,%g,[N],%g,[kgf]\n",Thrust,Thrust/9.8);
	fprintf(aerof_pt,"入力パワ,%g,[W]\n",P_input);
	fprintf(aerof_pt,"有効仕事,%g,[W]\n",Thrust * V_cruise);
	fprintf(aerof_pt,"プロペラ効率,%g,[-]\n",Eta);
	
	printf("\n%sに翼素データも書き込みますか[y/n]？	>",File.aero);
	y_n = getchar();
	gets(dummy);
	
	if(y_n == 'y'){
		fprintf(aerof_pt,"\n翼素データ\n");
		fprintf(aerof_pt,"y_p/R[-],F[-],Γ/(VR)[-],u/V[-],w/V[-],q/V[-],Re[-],phi[deg],α[deg],Cl[-],Cd[-],u - w * tanφ\n");
		for(j=2; j<=JX-1; j++){
			fprintf(aerof_pt,"%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n"
				,R_PT->y_p[j],ID_PT->loss_func[j],GAM[j],ID_PT->v.u[j],ID_PT->v.w[j],ID_PT->v.q[j],
				BS_DAT[j].Re_num,ID_PT->phi[j]*180/PI,Alpha_rad[j]*180/PI,
				BS_DAT[j].c.l,
				BS_DAT[j].c.d,
				ID_PT->v.u[j] - ID_PT->v.w[j]*tan(ID_PT->phi[j]));
		}
	}	
	fclose(aerof_pt);
}