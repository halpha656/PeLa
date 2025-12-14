/************************************************
		Spec_read_write.c
		
		 空力設計済みプロペラの
			翼型空力データファイル名、
			回転半径、設計内側半径、
			設計有効迎角、
			翼弦chord、捻り角分布β
		を読み書き
*************************************************/

#include"Prop_design.h"

//書き込む	Output()から呼ばれる　
void Write_spec(FILE *const aerof_pt, const struct profile_spec *const SPEC_PT)
{
	int j;

	//性能計算用データ書き込み
	fprintf(aerof_pt,"<性能計算・形状設計用データ>\n");
	fprintf(aerof_pt,"<翼型データ>\n");
	fprintf(aerof_pt,"%s\n",Profile_air.resource);

	fprintf(aerof_pt,"<スペック>\n");
	fprintf(aerof_pt,"進行率,%f\n",Adv);
	fprintf(aerof_pt,"設計有効迎角,%f\n",SPEC_PT->design_al_rad);
	fprintf(aerof_pt,"%f,%f\n",SPEC_PT->radius,SPEC_PT->r0);
	for(j=1; j<=JX-1; j++){
		fprintf(aerof_pt,"%f , %f\n",SPEC_PT->chord[j],SPEC_PT->beta[j]);
	}
	fprintf(aerof_pt,"<性能計算用データ終わり>");

}

//読み込む	自分で開ける閉じる
void Read_spec(struct profile_spec *const spec_PT,struct fixed_points *const r_PT)
{
	FILE *aerof_pt;
	char line[300];
	char resourcef[100];
	char *comma_pt;
	int j;
	
	do{
		printf("\n空力設計済みプロペラ名？\n>");
		scanf("%s",File.propname);	
	}while(Make_files() == FAILURE);
	gets(dummy);

	aerof_pt = fopen(File.aero,"r");
	while(fgets(line,sizeof(line),aerof_pt) != NULL){
		//データ開始を探す
		if(strstr(line,"<性能計算・形状設計用データ>") == NULL){
			continue;
		}
		else{//見つけたら
			
			//翼型データファイル名
			fgets(line,sizeof(line),aerof_pt);	//"<翼型データ>"読み捨て	
			fgets(line,sizeof(line),aerof_pt);
			sscanf(line,"%s",resourcef);	
			//フォーマットでとれないコンマ
			comma_pt = strchr(resourcef,',');
			if(comma_pt != NULL)	*comma_pt = '\0';		//あったらとる
				
			fgets(line,sizeof(line),aerof_pt);	//"<スペック>"読み捨て						
			//進行率
			fgets(line,sizeof(line),aerof_pt);
			sscanf(line,"進行率 , %lf",&Adv);

			//設計有効迎角
			fgets(line,sizeof(line),aerof_pt);
			sscanf(line,"設計有効迎角 , %lf",&spec_PT->design_al_rad);

			//回転半径、翼根半径
			fgets(line,sizeof(line),aerof_pt);	
			sscanf(line,"%lf , %lf",&(spec_PT->radius),&(spec_PT->r0));
			
			//翼弦、捩り角
			for(j=1; j<=JX-1; j++){
				fgets(line,sizeof(line),aerof_pt);
				sscanf(line,"%lf , %lf",&(spec_PT->chord[j]),&(spec_PT->beta[j]));
			}	
			
			//データ終了印確認
			fgets(line,sizeof(line),aerof_pt);
			if(strstr(line,"<性能計算用データ終わり>") != NULL){
				fclose(aerof_pt);
				
				//計算点設定
				Set_point(r_PT,spec_PT);
				//翼型データ読み込み
				Profile_air = Read_Airfoil_Data(PROFILE,resourcef,spec_PT);
				//終了
				return;
			}
			else{
				printf("データ数が合いません　分割数等を確認してください\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	//EOF
	printf("%sには読み込めるデータがありません\n",File.aero);
	exit(EXIT_FAILURE);		
}