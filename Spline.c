/*****************************************************************
			Spline.c
			非周期関数用３次スプライン補間
			服部さんの
*****************************************************************/
#include"Prop_design.h"

void make_table(const int N, const double x[], const double y[], double z[]);
void make_table2(const struct line refline, double table[]);

// x = t におけるyの補間値を返す
double Spline(const double t,		//x
			  const int N,			//参考点個数(x[],y[]のサイズ)
			  const double x[],		//参考点x
			  const double y[])		//参考点y
{
	int i,j,k;
	double d,h;
	double z[200];			//補間用テーブル
	double point;

	make_table(N,x,y,z);

	i=0;
	j=N-1;

	while(i<j){
		k=(i+j)/2;
		if( x[k] < t )	i=k+1;
		else	j=k;
	}

	if(i>0)	i--;
	h =x[i+1]-x[i];
	d =t-x[i];
	
	point = (((z[i+1]-z[i])*d/h+z[i]*3)*d+((y[i+1]-y[i])/h-(z[i]*2+z[i+1])*h))*d+y[i];

	return point;
}

//参考点のx、y座標から補間用テーブルz[]に値を代入　　理解不能
//x[0],y[0]が端点	Nは参考点の個数	
void make_table(const int N, const double x[], const double y[], double z[])
{
	int i;
	double t;
	static double h[200], d[200];	//大きめ

	if(N > 200){
		printf("Error!	 make_table()\n参考点の個数が多すぎます\n");
		exit(EXIT_FAILURE);
	}

	z[0] = 0;
	z[N-1] = 0;
	for(i=0; i<N-1; i++){
		h[i  ] =  x[i+1] - x[i];
		d[i+1] = (y[i+1] - y[i]) / h[i];
	}

	z[1]=d[2]-d[1]-h[0]*z[0];
	d[1]=(x[2]-x[0])*2;

	for(i=1; i<N-2; i++){
		t = h[i] / d[i];
		z[i+1] = d[i+2] - d[i+1] - z[i] * t;
		d[i+1] = 2 * (x[i+2] - x[i]) - h[i] * t;
	}
	
	z[N-2] -= h[N-2] * z[N-1];

	for(i=N-2; i>0; i--)
		z[i] = (z[i] - h[i]*z[i+1]) / d[i];

}

struct line Spline2(struct line spline, const struct line line_ref)
{
	int i,j,k,n;
	double d,h;
	double table[200];

	make_table2(line_ref, table);

	for(n=0; n<=spline.num_plot-1; n++){
		i=0;
		j=line_ref.num_plot-1;

		while(i<j){
			k=(i+j)/2;
			if( line_ref.x[k] < spline.x[n] )	i=k+1;
			else	j=k;
		}

		if(i>0)	i--;
		h =line_ref.x[i+1] - line_ref.x[i];
		d =spline.x[n] - line_ref.x[i];
		
		spline.y[n] = (((table[i+1]-table[i])*d/h+table[i]*3)*d+((line_ref.y[i+1]-line_ref.y[i])/h-(table[i]*2+table[i+1])*h))*d+line_ref.y[i];
	}
	return spline;
}

void make_table2(const struct line refline, double table[])
{
	int i;
	double t;
	static double h[200], d[200];	//大きめ

	if(refline.num_plot > 200){
		printf("Error!	 make_table2()\n参考点の個数が多すぎます\n");
		exit(EXIT_FAILURE);
	}

	table[0] = 0;
	table[refline.num_plot-1] = 0;
	for(i=0; i<refline.num_plot-1; i++){
		h[i  ] =  refline.x[i+1] - refline.x[i];
		d[i+1] = (refline.y[i+1] - refline.y[i]) / h[i];
	}

	table[1]=d[2]-d[1]-h[0]*table[0];
	d[1]=(refline.x[2]-refline.x[0])*2;

	for(i=1; i<refline.num_plot-2; i++){
		t = h[i] / d[i];
		table[i+1] = d[i+2] - d[i+1] - table[i] * t;
		d[i+1] = 2 * (refline.x[i+2] - refline.x[i]) - h[i] * t;
	}
	
	table[refline.num_plot-2] -= h[refline.num_plot-2] * table[refline.num_plot-1];

	for(i=refline.num_plot-2; i>0; i--)
		table[i] = (table[i] - h[i]*table[i+1]) / d[i];

}
