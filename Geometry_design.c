/*************************************************************
			Geometry_design.c
			rootのchordとβをスプライン補間
			桁に合わせて厚み増し
			リブ作り→textに各リブのデータを出力（CAD用）
*************************************************************/
#include"Prop_design.h"

//最終断面形状（有次元）
static struct rib_data{
	double r;				//半径位置[mm]
	double chord;			//翼弦長[mm]
	double beta;			//捩り角[deg]
	double thick_gain;		//オリジナルからの翼厚増加率[-]
	double shaft_r;			//桁半径[mm]
	double root_weight;		//root翼型混合率
};

static struct rib_data rib_root[RX];	//root	  0～RX-1		まで使う
static struct rib_data rib_prof[JX];	//profile Prof_top～JX-1 まで使う

int Input_parameter(struct profile_spec *const spec_PT);
void Profile_design(struct fixed_points *const r_PT, struct profile_spec *const spec_PT);
void Root_design(struct fixed_points *const r_PT, struct profile_spec *const spec_PT);
void Gain_thickness(void);
void Output_geometry_data(const struct fixed_points *const R_PT);
void Make_Script(void);

static int Shaft_x;				//x方向桁位置％
static int Prof_top;			//空力設計部分使い始め位置index
static double Plank;			//プランク厚み[mm]
static double Margin;			//プランク分からさらに必要な余裕(上下別)[mm]
static double Cut;				//前縁切り取り分[mm]
static double R1, D1, R2, D2;	//桁（回転半径[mm]、その位置の外径[mm]  2点）

void Geometry_design(struct fixed_points *const r_PT,struct profile_spec *const spec_PT)
{
	Profile_design(r_PT,spec_PT);
	
	while(Input_parameter(spec_PT) == FAILURE);
	Root_design(r_PT,spec_PT);
	Gain_thickness();
	
	Output_geometry_data(r_PT);
	Make_Script();
	printf("\n結果は%sに出力し,スクリプトファイル（AutoCAD LT2000用）%sを作りました\n",File.geometry,File.cad);
}

//空力設計部分
void Profile_design(struct fixed_points *const r_PT, struct profile_spec *const spec_PT)
{
	int i;

	//長さ単位->mm	角度->degree
	for(i=2; i<=JX-1; i++){
		rib_prof[i].r = r_PT->y_p[i] * spec_PT->radius * 1000;
		rib_prof[i].chord = spec_PT->chord[i] * spec_PT->radius * 1000;
		rib_prof[i].beta = spec_PT->beta[i] * 180 / PI;
		rib_prof[i].root_weight = 0.0;
	}
}

int Input_parameter(struct profile_spec *const spec_PT)
{	
	char root_resource[100];
	int y_n;

	printf("\n設計内側半径より内側のroot部分の形状設計をおこないます\n");
	if(calc_mode == ROOTDESIGN)	{
		printf("初めから:%d\nCanopus:%d\ntest:%d\n>",NEW,CANOPUS,TEST);
		scanf("%d",&design_mode);
		gets(dummy);
	}
	
	switch(design_mode){
	case NEW:
		printf("root最内半径[mm]？              >");	scanf("%lf",&rib_root[0].r);
		printf("そこでのchord[mm]？             >");	scanf("%lf",&rib_root[0].chord);
		printf("root翼型データファイル名？(.txt)>");	scanf("%s",root_resource);
		Root_air = Read_Airfoil_Data(ROOT,root_resource,spec_PT);
		if(strcmp(Root_air.resource,"unknown") == 0)	return FAILURE;

		printf("\n桁の半径位置[mm]と外径[mm]を2点指定してください\n");
		printf("1点目\n");
		printf("回転半径[mm]? >");		scanf("%lf",&R1);
		printf("桁外径[mm]?   >");		scanf("%lf",&D1);
		printf("2点目\n");
		printf("回転半径[mm]? >");		scanf("%lf",&R2);
		printf("桁外径[mm]?   >");		scanf("%lf",&D2);
		
		printf("\nプランク厚み[mm]?    >");	scanf("%lf",&Plank);
		printf("余裕[mm]?            >");	scanf("%lf",&Margin);
		printf("桁位置[%%]?           >");	scanf("%d",&Shaft_x);
		printf("前縁カット長さ[mm]?  >");	scanf("%lf",&Cut);
		gets(dummy);
		
		printf("入力ミス？	（ ある:y ない:n )	>");	y_n = getchar();	gets(dummy);
		if(y_n == 'y')	return FAILURE;
		else	break;

	case CANOPUS:
		rib_root[0].r = 50.0;
		rib_root[0].chord = 80.0;
		//桁はr=1500から5mmはみだす		tip側からとる
		R1 = 50.0;		D1 = 24.10;
		R2 = 1505.0;	D2 = 6.15;
		Plank = 1.5;
		Margin = 0.7;
		Shaft_x = 36;
		Cut = 9.0;
		Root_air = Read_Airfoil_Data(ROOT,"e856.txt",spec_PT);
		if(strcmp(Root_air.resource,"unknown") == 0){
			printf("Error!\nGeometry_design.c\nInput_parameter\ncase TEST\n");
			exit(EXIT_FAILURE);
		}
		else	break;

	case TEST:
		rib_root[0].r = 50.0;
		rib_root[0].chord = 80.0;
		R1 = 10.0;		D1 = 25.0;
		R2 = 1000.0;	D2 = 10.0;
		Plank = 1.0;
		Margin = 1.0;
		Shaft_x = 36;
		Cut = 10.0;
		Root_air = Read_Airfoil_Data(ROOT,"dae41.txt",spec_PT);
		if(strcmp(Root_air.resource,"unknown") == 0){
			printf("Error!\nGeometry_design.c\nInput_parameter\ncase TEST\n");
			exit(EXIT_FAILURE);
		}
		else	break;

	default:
		printf("Error!\nGeometry_design.c\nInput_parameter\ndesign_mode = %d\n",design_mode);
		exit(EXIT_FAILURE);
	}
	
	return 0;

}

//スプライン補間でrootの翼弦長と捩り角を決める
void Root_design(struct fixed_points *const r_PT, struct profile_spec *const spec_PT)
{
	//補間用
	double xx[JX+1];
	double yy1[JX+1], yy2[JX+1];
	
	double dr;		//root分割幅[mm]
	int i;
	
	rib_root[0].beta = (atan(Adv * 1000 / rib_root[0].r) + spec_PT->design_al_rad) * 180.0 / PI;
	
	//端点
	xx[0] = rib_root[0].r;
	yy1[0] = rib_root[0].chord;
	yy2[0] = rib_root[0].beta;
	
	//空力設計部分を全部使うとくびれるので、Prof_top～JX-1を使う
	for(Prof_top=2; Prof_top<=JX-1; Prof_top++){
		printf("Prof_top = %02d\n",Prof_top);		
		dr = ( rib_prof[Prof_top].r - rib_root[0].r ) / RX;

		for(i=1; i<=JX-Prof_top; i++){
			xx[i] = rib_prof[i-1+Prof_top].r;
			yy1[i] = rib_prof[i-1+Prof_top].chord;
			yy2[i] = rib_prof[i-1+Prof_top].beta;
		}

		//翼弦長、捩り角を補間
		for(i=1; i<=RX-1; i++){
			rib_root[i].r = rib_root[0].r + i * dr;
			rib_root[i].chord = Spline(rib_root[i].r,JX+1-Prof_top,xx,yy1);
			rib_root[i].beta  = Spline(rib_root[i].r,JX+1-Prof_top,xx,yy2);
			printf("chord[%02d] = %g\tbeta[%02d] = %g\n",i,rib_root[i].chord,i,rib_root[i].beta);
		}
		//くびれてなかったら出る
		if(rib_root[0].chord < rib_root[1].chord)	break;
	}

	//chordからroot側翼型の割合を決める
	for(i=0; i<=RX-1; i++){
		rib_root[i].root_weight = Linear_Polation(rib_root[i].chord, rib_root[0].chord, rib_prof[Prof_top].chord, 1.0, 0.0);
	}

	printf("補間で無視した空力設計部分 -> 半径%g%%～%g%%\n",r_PT->y_p[2]*100,r_PT->y_p[Prof_top-1]*100);

}

//桁がはみ出ないようにリブを厚み増し
void Gain_thickness(void)
{
	double t_shaft;			//桁位置での厚み[mm]
	double root_t_ratio;	//root側　オリジナル翼型の厚み比[-](桁位置）
	double prof_t_ratio;	//profile側
	int i;

	root_t_ratio = Root_air.u.y[Shaft_x]    - Root_air.l.y[Shaft_x];
	prof_t_ratio = Profile_air.u.y[Shaft_x] - Profile_air.l.y[Shaft_x];

	//root
	for(i=0; i<=RX-1; i++){
		rib_root[i].shaft_r = Linear_Polation(rib_root[i].r,R1,R2,D1,D2) / 2.0;

		t_shaft = (rib_root[i].root_weight * root_t_ratio + (1.0 - rib_root[i].root_weight) * prof_t_ratio) * rib_root[i].chord;
		
		if(t_shaft/2.0 > rib_root[i].shaft_r + Plank + Margin)	rib_root[i].thick_gain = 1.0;
		else	rib_root[i].thick_gain = ( rib_root[i].shaft_r + Plank + Margin ) * 2.0 / t_shaft;
	}
	//profile
	for(i=Prof_top; i<=JX-1; i++){
		rib_prof[i].shaft_r = Linear_Polation(rib_prof[i].r,R1,R2,D1,D2) / 2.0;
		if(rib_prof[i].shaft_r <= 0.0)	rib_prof[i].shaft_r = 0.1;
		t_shaft = prof_t_ratio * rib_prof[i].chord;
		
		if(t_shaft/2.0 > rib_prof[i].shaft_r + Plank + Margin)	rib_prof[i].thick_gain = 1.0;
		else	rib_prof[i].thick_gain = ( rib_prof[i].shaft_r + Plank + Margin ) * 2.0 / t_shaft;
	}
}		

//形状データ出力
void Output_geometry_data(const struct fixed_points *const R_PT)
{
	FILE *geometryf_pt;
	int i;

	geometryf_pt = fopen(File.geometry,"w");

	fprintf(geometryf_pt,"\n\n>>>>>>>>>>>>>>>>>>>>>>rootつき形状データ>>>>>>>>>>>>>>>>>>>\n");
	fprintf(geometryf_pt,"Root側翼型,%s\n",Root_air.name);
	fprintf(geometryf_pt,"補間で無視した空力設計部分 -> 半径%g%%～%g%%\n",R_PT->y_p[2]*100,R_PT->y_p[Prof_top-1]*100);
	fprintf(geometryf_pt,"プランク厚み,%g,[mm]\n",Plank);
	fprintf(geometryf_pt,"余裕,%g,[mm]\n",Margin);
	fprintf(geometryf_pt,"前縁カット,%g,[mm]\n",Cut);
	fprintf(geometryf_pt,"桁位置,%d,[%%]\n",Shaft_x);
	fprintf(geometryf_pt,"桁形状,(r=%g、φ%g),(r=%g、φ%g)\n",R1,D1,R2,D2);
	fprintf(geometryf_pt,",リブ番号,半径[mm],翼弦長[mm],捩り角[deg],桁外半径[mm],翼厚増分[-],root側翼型の割合[-]\n");
	fprintf(geometryf_pt,"<Root>");
	for(i=0; i<=RX-1; i++){
		fprintf(geometryf_pt,",R%02d,%g,%g,%g,%g,%g,%g\n",
				i,rib_root[i].r,rib_root[i].chord,rib_root[i].beta,
				rib_root[i].shaft_r,rib_root[i].thick_gain,rib_root[i].root_weight);
	}
	fprintf(geometryf_pt,"<Profile>");
	for(i=Prof_top; i<=JX-1; i++){
		fprintf(geometryf_pt,",P%02d,%g,%g,%g,%g,%g,%g\n",
				i,rib_prof[i].r,rib_prof[i].chord,rib_prof[i].beta,
				rib_prof[i].shaft_r,rib_prof[i].thick_gain,rib_prof[i].root_weight);
	}

	fclose(geometryf_pt);
}

//Auto CAD LT 2000 対応スクリプトファイル作成
void Make_Script(void)
{
	FILE *scrf_pt;
	double root_y_shaft,prof_y_shaft;
	double c, w;
	double x_shaft, y_shaft;
	int i,j;
	
	scrf_pt = fopen(File.cad,"w");

	root_y_shaft = ( Root_air.u.y[Shaft_x]    + Root_air.l.y[Shaft_x]    ) / 2.0;
	prof_y_shaft = ( Profile_air.u.y[Shaft_x] + Profile_air.l.y[Shaft_x] ) / 2.0;
	
	//オブジェクトスナップ//**できれば自動でオフ
	fprintf(scrf_pt,"dsettings\n");

	fprintf(scrf_pt,"view swiso\n");
	fprintf(scrf_pt,"layer new sub color yellow sub lw 0.0 sub set sub\n\n");
	fprintf(scrf_pt,"xline h 0,0,0\n\n");
	fprintf(scrf_pt,"xline v 0,0,0\n\n");
	
	//Root
	fprintf(scrf_pt,"layer new R-1 set R-1 freeze sub\n\n");	//dummy layer
	for(i=0; i<=RX-1; i++){
		c = rib_root[i].chord;
		w = rib_root[i].root_weight;
		//桁中心の座標
		x_shaft = (double)Shaft_x / 100 * c;
		y_shaft = (root_y_shaft * w + prof_y_shaft * ( 1.0 - w ) )
					* c * rib_root[i].thick_gain;
		
		//layer作成
		fprintf(scrf_pt,"layer new R%02d color blue R%02d lw 0.0 R%02d set R%02d freeze R%02d \n",i,i,i,i,i-1);
		//断面形
		fprintf(scrf_pt,"spline\n");
			for(j=100;j>=0;j--){
				fprintf(scrf_pt,"%g,%g,0\n"
					,Root_air.u.x[j] * c
					,(Root_air.u.y[j] * w + Profile_air.u.y[j] * (1.0 - w)) * c * rib_root[i].thick_gain);
			}
			for(j=0; j<=100; j++){
				fprintf(scrf_pt,"%g,%g,0\n"
					,Root_air.l.x[j] * c
					,(Root_air.l.y[j] * w + Profile_air.l.y[j] * (1.0 - w)) * c * rib_root[i].thick_gain);
			}
			fprintf(scrf_pt,"\n0 0\n");
			//プランク分内側へオフセット
			fprintf(scrf_pt,"offset %g 0,0,0\n%g,%g,0\n\n",Plank,x_shaft,y_shaft);
		//前縁切り取り線
		fprintf(scrf_pt,"xline v %g,0,0 \n",Cut);
		//桁外型
		fprintf(scrf_pt,"circle %g,%g,0 %g\n",x_shaft,y_shaft,rib_root[i].shaft_r);
		//翼弦
		fprintf(scrf_pt,"color red\n");
		fprintf(scrf_pt,"xline h 0,0,0\n\n");
		//前縁,後縁
		fprintf(scrf_pt,"xline v 0,0,0 %g,0,0 \n",rib_root[i].chord);
		//リブ番号
		fprintf(scrf_pt,"color magenta\n");
		fprintf(scrf_pt,"text j ml %g,%g,0 4.5 0 R%02d\n",0.55*rib_root[i].chord,y_shaft,i);
		//その他
		fprintf(scrf_pt,"color green\n");
		fprintf(scrf_pt,"circle %g,%g,0 %g\n",x_shaft,y_shaft,rib_root[i].shaft_r + 0.5);
		fprintf(scrf_pt,"xline h %g,%g,0 \n",x_shaft,y_shaft);
		
		//移動（原点：前縁→桁中心）,回転
		fprintf(scrf_pt,"ai_selall move %g,%g,0 0,0,0\n",x_shaft,y_shaft);
		fprintf(scrf_pt,"ai_selall rotate 0,0,0 %g\n",-rib_root[i].beta);
		//3D
		fprintf(scrf_pt,"ai_selall move 0,0,0 0,0,%g\n",-rib_root[i].r);
		fprintf(scrf_pt,"color bylayer\n");
	}
	
	//Proflile
	fprintf(scrf_pt,"layer new P%02d set P%02d \n",Prof_top-1,Prof_top-1);	//dummy
	fprintf(scrf_pt,"layer freeze R%02d \n",RX-1);
	for(i=Prof_top; i<=JX-1; i++){
		c = rib_prof[i].chord;
		//桁中心の座標
		x_shaft = (double)Shaft_x / 100 * c;
		y_shaft = prof_y_shaft * c * rib_prof[i].thick_gain;
		
		fprintf(scrf_pt,"layer new P%02d color blue P%02d lw 0.0 P%02d set P%02d freeze P%02d \n",i,i,i,i,i-1);
		
		fprintf(scrf_pt,"spline\n");
			//upper
			for(j=100; j>=0; j--){
				fprintf(scrf_pt,"%g,%g,0\n",Profile_air.u.x[j] * c,Profile_air.u.y[j] * c * rib_prof[i].thick_gain);
			}
			//lower
			for(j=0; j<=100; j++){
				fprintf(scrf_pt,"%g,%g,0\n",Profile_air.l.x[j] * c,Profile_air.l.y[j] * c * rib_prof[i].thick_gain);
			}
			fprintf(scrf_pt,"\n0 0\n");
			fprintf(scrf_pt,"offset %g 0,0,0\n%g,%g,0\n\n",Plank,x_shaft,y_shaft);
		fprintf(scrf_pt,"xline v %g,0,0 \n",Cut);
		fprintf(scrf_pt,"circle %g,%g,0 %g\n",x_shaft,y_shaft,rib_prof[i].shaft_r);
		
		fprintf(scrf_pt,"color red\n");
		fprintf(scrf_pt,"xline h 0,0,0\n\n");
		fprintf(scrf_pt,"xline v 0,0,0 %g,0,0 \n",rib_prof[i].chord);
		
		fprintf(scrf_pt,"color magenta\n");
		fprintf(scrf_pt,"text j ml %g,%g,0 4.5 0 P%02d\n",0.55*rib_prof[i].chord,y_shaft,i);

		fprintf(scrf_pt,"color green\n");
		fprintf(scrf_pt,"circle %g,%g,0 %g\n",x_shaft,y_shaft,rib_prof[i].shaft_r + 0.5);
		fprintf(scrf_pt,"xline h %g,%g,0 \n",x_shaft,y_shaft);		
		
		fprintf(scrf_pt,"ai_selall move %g,%g,0 0,0,0\n",x_shaft,y_shaft);
		fprintf(scrf_pt,"ai_selall rotate 0,0,0 %g\n",-rib_prof[i].beta);
		fprintf(scrf_pt,"ai_selall move 0,0,0 0,0,%g\n",-rib_prof[i].r);
		fprintf(scrf_pt,"color bylayer\n");
	}

	//全画層freeze解凍
	fprintf(scrf_pt,"layer thaw ");
		for(i=0; i<=RX-1; i++)fprintf(scrf_pt,"R%02d,",i);
		for(i=Prof_top; i<=JX-1; i++)fprintf(scrf_pt,"P%02d,",i);
		fprintf(scrf_pt,"sub,R-1,P%02d\n\n",Prof_top-1);
	
	//上向きにする
	fprintf(scrf_pt,"ucs new za 0,0,0 0,1,0\n");
	fprintf(scrf_pt,"ai_selall rotate 0,0,0 180\n");
	fprintf(scrf_pt,"ucs world\n");

	//オブジェクトスナップ
	fprintf(scrf_pt,"dsettings\n");
	
	fclose(scrf_pt);

	
}