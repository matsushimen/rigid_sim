//
//  mol_def.h
//  
//
//  Created by 松島　佑樹 on 2016/06/16.
//
//

#ifndef _mol_def_h
#define _mol_def_h

#define DEMENSION 3
#define PI (M_PI)
#define T (50.0)
#define k_B (1.0)//ボルツマン定数
#define N_A (6.02214129e+23)//アボガドロ定数
#define ANG (1.0e-10)
#define NENSEI (1.0/6.0)


//座標の数値化
typedef enum type_coordinate
{
    X, Y, Z
} COORDINATE;


//粒子の構造体
typedef struct molecule
{
	unsigned int number;//index of molecule
    double r;//radius of molecule
	double p_0[DEMENSION];//first position of molecule
	double p[DEMENSION];//position of molecule
	double p_n[DEMENSION];//new position of molecule
	double v[DEMENSION];//velocity of molecule
	double v_n[DEMENSION];//new velocity of molecule
    double m;//mass of molecule
    double f[DEMENSION];
    double gamma;
    double scale;
    double def_dist_a;
    double def_dist_b;
    int bond_pair;
    int breakpair1_b;
    int breakpair1_a;
    int breakpair2_b;
    int breakpair2_a;
} Molecule;

Molecule *mlcl = NULL;

#endif
