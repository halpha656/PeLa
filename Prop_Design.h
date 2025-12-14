/*******************************************************
			Prop_design.h
*******************************************************/
#define PI   3.141592654
#define JX   51			//ブレード分割数
#define J70  30			//だいたい７０％位置
#define RX   20			//root分割数
#define IX   1300		//後流渦分割数
#define DX_0 1.0e-6
#define SX   1.01
#define CONVERGE_RLX	0.001
#define CONVERGE_LAG	0.01
#define FAILURE 1

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

//計算用の固定点座標（プロペラの回転半径で無次元化）
struct fixed_points{
	double y[JX+1];			//ブレード上積分点y座標	
	double dy[JX+1];		//分割幅（yの変化量）
	double y_p[JX+1];		//control point　y座標
	double x_v[IX+1];		//後流渦の計算点x座標
	double accel_mdl[IX+1];	//後流渦の加速っぷりを近似		//**人力なら要らない
};

//[j][k]-->j番目の循環がk番目のcontrol pointに与える影響
struct factor{
	double a0[JX+1][JX+1];	//z方向（回転方向)	//着目してるブレード自身から
	double a1[JX+1][JX+1];						//もう一本のブレードから
	double b0[JX+1][JX+1];	//x方向（機体進行方向）
	double b1[JX+1][JX+1];
};

//局所流入速度（機体の巡航速度で無次元化）
struct velocity{
	double u[JX+1];			//x方向誘導速度
	double w[JX+1];			//z方向誘導速度
	double q[JX+1];			//局所流入速度
};

struct induced{
	struct factor fc;
	struct velocity v;		
	double phi[JX+1];		//流入角[rad]
	double loss_func[JX+1];	//速度減少関数(A&Lの性能計算でのみ使用)
	double disp_v_half;		//渦の移動速度（ブレード上）
};
//空力設計上の形
struct profile_spec{
	double radius;			//回転半径[m]
	double r0;				//翼根半径[m]
	double design_al_rad;	//設計有効迎角[rad]
	double chord[JX+1];		//翼弦長[-]（回転半径で無次元化）
	double beta[JX+1];		//捩り角[rad]
};

//プロペラ全体の抵抗係数（負の推力）とトルク係数
struct Coefs{
	double dg[3];
	double dgtotal;		//dg[0]+dg[1]+dg[2]
	double tq[3];
	double tqtotal;		//tq[0]+tq[1]+tq[2]
};

//局所的な空力係数
struct blade_section_coefs{
	double l;				//Cl
	double l_opt;			//設計Cl（補間値）
	double dm[3];			//Cd=dm[0]+dm[1]*Cl+dm[2]*Cl^2
	double d;				//Cd
};

//翼素のデータ
struct blade_section_data{
	struct blade_section_coefs c;	//空力係数
	double Re_num;					//レイノルズ数
	double stall;					//失速角[rad]
};

//polar, Cl-α
struct character_crvs{
	double Re;				//ポーラのレイノルズ数
	double cl_design;		//設計Cl　迎角からReに対して決める
	double cd[20];			//Cl=index/10.0に対するCd (polar)
	double grad_cl_al;		//傾き 1/[rad]
	double cl_al0;			//cl-α切片[-]
	double stall_rad;		//失速角[rad]
};

//曲線
struct line{
	double x[101];		//x座標
	double y[101];		//y座標
	int num_plot;		//点の数
};

//翼型データ textファイルから読み取る
struct airfoil_data{
	char resource[20];				//データファイル名(.txt)
	char name[20];					//翼型名
	int num_crvs;					//Reのパターン数
	int polar_size;					//polarの最大Cl
	struct character_crvs crvs[10];	//polarとCl-α(Re別)
	struct line u;					//upper
	struct line l;					//lower
};

struct range{
	double lower;	//下限
	double upper;	//上限
};

//ループの回数をカウント,時間を記録
struct loop_counter{
	int relax;
	int lag;
	double time_rec;
};

struct resultf{
	char propname[30];	//プロペラの名前
	char aero[40];		//空力設計部データ・性能計算結果出力
	char geometry[40];	//形状データ出力
	char cad[40];		//AutoCADのスクリプト
};

//循環分布を設定・更新するときのモード
enum gam_change_mode{INITIALIZE,	//初期化
					 RELAXATION,	//偏微分を用いた最適化（トルクのウエイト一定のもとで）
					 LAGLANGE}		//全体を一定の割合(κ倍)で変えてトルクを設計値に合わせる
					gcmode;

//粘性（形状抵抗）の有無{非粘性, 粘性有り}
enum viscous_contribution{INVISCID, VISCID}visc;

//{設計, root設計, 設計点外性能計算}
enum caliculation_mode{DESIGN, ROOTDESIGN, OFFDESIGN}calc_mode;

enum design_mode{NEW, CANOPUS, TEST}design_mode;

enum airfoil{PROFILE, ROOT};

//結果出力ファイル
//FILE *Resultf_pt;
//extern char Resultf[30];
extern struct resultf File;				//結果出力ファイル群

extern struct loop_counter total_lp;
extern struct loop_counter invis_lp;

extern double Re_prop;				//基準レイノルズ数	R*V/ν[-]
extern double Adv;					//進行率	V/(R*Ω)[-]
extern double CT_target;			//設定トルク係数[-]
extern double Alpha_rad[JX+1];		//有効迎角[rad]
extern double Pitch_rad;			//ピッチ角[rad]

extern struct airfoil_data Profile_air;	//空力設計翼型データ
extern struct airfoil_data Root_air;	//Root翼型データ

extern char dummy[10];				//改行コード除去

extern void Set_point(struct fixed_points *const r_PT, const struct profile_spec *const SPEC_PT);
extern void Get_blade_section_data(struct blade_section_data bs_dat[],const double GAM[], const struct induced *const ID_PT,const struct profile_spec *const SPEC_PT);
extern void Change_Gammas(enum gam_change_mode gcmode,const struct fixed_points *const R_PT,const double d_gam[JX+1],const double KAPPA,double gamma[],struct induced *const id_PT,struct blade_section_data bs_dat[],struct profile_spec *const spec_PT);
extern void Evaluate(const struct fixed_points *const R_PT,const struct induced *const ID_PT,const double GAM[],const struct blade_section_data BS_DAT[],struct profile_spec *const spec_PT);
extern void Biot_Savart(struct induced *const ID_pt,const struct fixed_points *const R_pt,const double GAM[]);
extern struct Coefs Coefficients(const struct blade_section_data BS_DAT[],const double GAM[],const struct fixed_points *const R_PT,const struct induced *const ID_PT,const struct profile_spec *const SPEC_PT);
extern struct airfoil_data Read_Airfoil_Data(enum airfoil read_mode, const char resourcef[],const struct profile_spec *const SPEC_PT);
extern void Read_spec(struct profile_spec *const spec_PT,struct fixed_points *const r_PT);
extern double Spline(const double t,const int N,const double x[],const double y[]);
extern struct line Spline2(struct line spline, const struct line line_ref);
extern double Linear_Polation(const double x, const double xref1, const double xref2,const double yref1, const double yref2);
extern int Make_files(void);