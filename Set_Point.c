/********************************************************
		Set_Point.c
		固定点の座標を設定
*********************************************************/

#include"Prop_design.h"

//計算点の座標を設定 以後変更無し
void Set_point(struct fixed_points *const r_PT, const struct profile_spec *const SPEC_PT)
{
	int i, j, k;
	double y0;
	double dx_v = DX_0;
	double theta[JX+1];

	y0 = SPEC_PT->r0 / SPEC_PT->radius;

	r_PT->x_v[0] = 0.0;
	r_PT->accel_mdl[0] = 1.0;

	//vortex sheet x座標設定
	//**渦の移動速度の増加(1 + disp_v_half --> 1 + 2*disp_v_half)を近似する関数を作る	人力ならいらない
	for(i=1; i<=IX; i++){
		dx_v = dx_v * SX;
		r_PT->x_v[i] = r_PT->x_v[i-1] + dx_v;
		r_PT->accel_mdl[i] = 1.0;		//2.0 - exp( - r_PT->x_v[i]/2.0 );
	}

	//control point 設定（cosine分割）
	for(j=1; j<=JX; j++){
		theta[j] = (double)( j - 1 ) / (JX - 1 ) * PI ;		
		r_PT->y_p[j] = y0 + ( 1.0 - y0 ) / 2 - ( 1.0 - y0 ) / 2 * cos(theta[j]);
	}

	//integration point 設定
	for(k=1; k<=JX-1; k++){
		r_PT->y[k] = y0 + ( 1.0 - y0 ) / 2 - ( 1.0 - y0 ) / 2 * cos(theta[k] + PI / (( JX - 1 )*2.0) );
	}

	//dy[]
	for(k=2; k<=JX-1; k++){
		r_PT->dy[k] = r_PT->y[k] - r_PT->y[k-1];
	}
	r_PT->dy[1] = r_PT->dy[2];

}