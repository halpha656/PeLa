/*************************************************************
		規格化済みの翼型データを読み込む
*************************************************************/
//読み込み可能ファイルが指定された場合、air.resource に resourcef をコピー
//読み込めなかった場合、air.resource に "unknown"　を書く
//**エラーチェックがややこしい
#include"Prop_design.h"

struct airfoil_data Read_Airfoil_Data(enum airfoil read_mode, const char resourcef[],const struct profile_spec *const SPEC_PT)
{
	FILE *resourcef_pt;
	char fname[100];
	char line[300];
	int n = 0;
	int min = 100;
	int data_num;
	struct airfoil_data air;
	struct line line_u0, line_u, line_l;
	int i;
	
	for(i=0; i<=100; i++){
		air.u.x[i] = (double)i / 100.0;
		air.l.x[i] = (double)i / 100.0;
	}
	
	switch(read_mode){
	case PROFILE:
		sprintf(fname,".\\airfoil\\%s",resourcef);
		break;
	case ROOT:
		sprintf(fname,".\\airfoil\\root\\%s",resourcef);
		break;
	}

	resourcef_pt = fopen(fname,"r");
	if(resourcef_pt == NULL){
		printf("%sを開くことができません。\n",fname);
		strcpy(air.resource,"unknown");
		return air;
	}
	else{
		//規格化印を確認
		fgets(line,sizeof(line),resourcef_pt);
		if(strncmp(line,"standardized for propeller design",33) != 0){
			printf("%sは規格化されていません。データを読み込めません。\n",fname);
			fclose(resourcef_pt);
			strcpy(air.resource,"unknown");
			return air;
		}
	}

	//読み込み可能ファイルが指定された場合
	strcpy(air.resource,resourcef);
	while(fgets(line,sizeof(line),resourcef_pt) != NULL){

		//１列目に'>'があったらデータを読み込む　それ以外の行は無視
		if(line[0] == '>'){
			//翼型名">name"
			if(strncmp(line,">name",5) == 0 )	sscanf(line,">name %s",&air.name);
			//レイノルズ数">Re"
			else if( strncmp(line,">Re",3) == 0 ){
				++n;
				sscanf(line,">Re %lf",&air.crvs[n].Re);
			}
			//配列番号/10に対応するCd">Cd"
			else if( strncmp(line,">Cd",3) == 0 ){

				data_num = sscanf(line,">Cd %lf ,  %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf ,",
								&air.crvs[n].cd[0],&air.crvs[n].cd[1],&air.crvs[n].cd[2],&air.crvs[n].cd[3],&air.crvs[n].cd[4],&air.crvs[n].cd[5],&air.crvs[n].cd[6],
								&air.crvs[n].cd[7],&air.crvs[n].cd[8],&air.crvs[n].cd[9],&air.crvs[n].cd[10],&air.crvs[n].cd[11],&air.crvs[n].cd[12],&air.crvs[n].cd[13],
								&air.crvs[n].cd[14],&air.crvs[n].cd[15]);
				if(data_num < min )	min = data_num;
			}
			//Cl-αの傾き[1/rad]">grad"
			else if( strncmp(line,">grad",5) == 0 )	sscanf(line,">grad %lf",&air.crvs[n].grad_cl_al);
			//Cl-αの切片">Cl0"
			else if( strncmp(line,">Cl0",4) == 0 )	sscanf(line,">Cl0 %lf",&air.crvs[n].cl_al0);
			//失速迎角[rad]">stall"
			else if(strncmp(line,">stall",6) == 0 )	sscanf(line,">stall %lf",&air.crvs[n].stall_rad);
			
			//形状データ（座標 スペース区切り
			else if(strncmp(line,">plot",5) == 0 ){
				//上面座標読み込み	（後縁→前縁）
				fgets(line,sizeof(line),resourcef_pt);	sscanf(line,"%lf %lf",&line_u0.x[0],&line_u0.y[0]);
				fgets(line,sizeof(line),resourcef_pt);	sscanf(line,"%lf %lf",&line_u0.x[1],&line_u0.y[1]);
				for(i=1; line_u0.x[i] < line_u0.x[i-1]; i++){
					fgets(line,sizeof(line),resourcef_pt);
					sscanf(line,"%lf %lf",&line_u0.x[i+1],&line_u0.y[i+1]);
				}
				line_l.x[0] = 0.0;	line_l.y[0] = 0.0;
				line_l.x[1] = line_u0.x[i];	line_l.y[1] = line_u0.y[i];
				
				if(line_u0.x[i-1] == 0.0)	line_u0.num_plot = i;
				else{
					line_u0.x[i] = 0.0;		line_u0.num_plot = i+1;
				}

				//下面（前縁→後縁）
				fgets(line,sizeof(line),resourcef_pt);
				for(i=2; sscanf(line,"%lf %lf",&line_l.x[i],&line_l.y[i])==2; i++){
					fgets(line,sizeof(line),resourcef_pt);
				}
				line_l.num_plot = i;
				//上面はスプライン用にひっくり返す（前縁→後縁）
				line_u.num_plot = line_u0.num_plot;	
				for(i=0; i<=line_u.num_plot-1; i++){
					line_u.x[i] = line_u0.x[line_u0.num_plot-1-i];
					line_u.y[i] = line_u0.y[line_u0.num_plot-1-i];
				}
			}
				
		}
	}
	fclose(resourcef_pt);
	
	air.num_crvs = n;
	air.polar_size = min - 1;

	//スプライン補間でx方向に100分割した座標になおす
	air.u.num_plot = 101;	air.u = Spline2(air.u, line_u);
	air.l.num_plot = 101;	air.l = Spline2(air.l, line_l);
	
	for(n=1; n<=air.num_crvs; n++){
		air.crvs[n].cl_design = air.crvs[n].grad_cl_al * SPEC_PT->design_al_rad + air.crvs[n].cl_al0;
	}

	return air;

}
