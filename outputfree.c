//
//  rigid
//
//
//  Created by 松島　佑樹 on 2016/05/09.
//
//

#include <stdlib.h>
#include <stdio.h>
#include "mersenneTwister.h"
#include "myfunc.h"
#include <time.h>
#include <math.h>
#include <string.h>
#include <omp.h>


#define STEP (10000000000)
#define D (2.0e+5)//刻み数
#define EXVE (10000.0)//排除体積バネ定数
#define LINEAR (10000.0)//線形バネ接続
#define RIGID (10000.0)

////////座標の計算に必要な変数
unsigned int Number;//number of molecule
double  Max = 0.1,min = 100000000.0, sum = 0.0;
int i, j, n, m, k, read = 0, save = 0, large, small;
unsigned long long step;
double t = 1.0/(float)D, d, d_0;
FILE *gp;
time_t times;

//////ファイルを読み込む関数
void readposition(int argc, char *argv[]){
    
	if(read == 0){
		if (argc != 2){
			printf("\n\n    error : input position file\n\n");
			exit( 1 );
		}
        
		unsigned int i, j;
		char filename[128];
		FILE *fp;
        
		////ファイルを開く
		sprintf(filename, "%s", argv[1]);
        
		//ファイルが開けなかったら
		if ((fp = fopen(filename, "r")) == NULL){
            printf("\n\n    error : can not open %s to read\n\n", filename);
            exit( 1 );
		}
        
		//read molecule numbers
		if(fscanf(fp, "%u\n", &Number)!=1){
            exit( 1 );
        }
        
		//メモリを確保する
		mlcl = (Molecule *)malloc(Number * sizeof(Molecule));
        
		//メモリが確保できなかったら
		if(mlcl == NULL){
            printf("\n\n    error : can not sequre the memories of molecule\n\n");
            exit( 1 );
		}
        //大きい粒子の数を読み取る
        
        if(fscanf(fp, "%u\n", &large)!=1){
            exit( 1 );
        }
        //小さい粒子の数を読み取る
        if(fscanf(fp, "%u\n", &small)!=1){
            exit( 1 );
        }
        //最後にシミュレーションしたときのステップ
        if(fscanf(fp, "%lld\n", &step)!=1){
            exit(1);
        }
        //最後の速度を使う
        if(fscanf(fp, "%d\n", &save)!=1){
            exit(1);
        }
		//粒子の初期位置と半径を読み取る
		for(i = 0; i < Number; i++){
            
            if(fscanf(fp, "%u %lf %lf %lf %lf %lf %lf\n",&mlcl[i].number, &mlcl[i].p[X], &mlcl[i].p[Y], &mlcl[i].p[Z], &mlcl[i].r, &mlcl[i].m, &mlcl[i].scale)!=7){
                exit( 1 );
            }
            mlcl[i].gamma = 1.0;
            mlcl[i].breakpair1_a = (i + 1) % Number;
            mlcl[i].breakpair1_b = (Number + i - 1) % Number;
            mlcl[i].breakpair2_a = (i + 2) % Number;
            mlcl[i].breakpair2_b = (Number + i -2) % Number;
            mlcl[i].bond_pair = -1;
            mlcl[i].f[X] = 0.0;
            mlcl[i].f[Y] = 0.0;
            mlcl[i].f[Z] = 0.0;

		}
        
		////ファイルを閉じる
		fclose(fp);
		read = 1;
        printf("\n");
	}
}

void molsposition(void){
    
	double p1,p2,theta,psi;
	double Fx[Number],Fy[Number],Fz[Number];
    double count[Number][Number];
    double list_no[Number], list_pair[Number][Number];
    double r_rist = 4.0;
    double d;
    
   

    
	for(step = step; step < STEP; step++){
        if(save == 0){//初期位置の保存と初速度の付加
            
            for(i = 0; i < Number; i++){
                mlcl[i].p_n[X] = mlcl[i].p[X];
                mlcl[i].p_n[Y] = mlcl[i].p[Y];
                mlcl[i].p_n[Z] = mlcl[i].p[Z];
                mlcl[i].f[X] = 0.0;
                mlcl[i].f[Y] = 0.0;
                mlcl[i].f[Z] = 0.0;
                
            }
            for(i = 0; i < Number; i++){//count配列の初期化
                for(j = 0; j < Number; j++)
                {
                    count[i][j] = 0.0;
                }
            }
            save = 1;
        }
        if(save == 2)
        {
            FILE *pp;
            char matrixfile[128];
            sprintf(matrixfile, "resultfree/breakcountmatrix%lld.txt",step);
            if ((pp = fopen(matrixfile, "r")) == NULL)
            {
                printf("\n\n    error : can not open %s to read\n\n", matrixfile);
                exit( 1 );
            }
            for(i=0;i<Number;i++){
                
                for(j=0;j<Number;j++){
                    
                    if ((fscanf(pp,"%lf",&count[i][j])) != 1){
                        
                        printf("\n\n    error : can not read %s to read\n\n", matrixfile);
                        exit( 1 );
                    }
                    count[i][j] = count[i][j]*step;
                }
            }
            save = 1;
            fclose(pp);
        }
		else
		{
            /////////リスト更新///////////
            if(step % 50 == 0)//50ステップごとにリストを作成
            {
                for(k = 0; k < Number; k++)
                {
                    m = 0;
                    for(n = 0; n < Number; n++)
                    {
                        d = dist(k,n);
                        if((d < r_rist) && (k!=n))
                        {
                            m++;
                            list_no[k] = m;
                            list_pair[k][m] = n;
                        }
                    }
                    if(m == 0)
                    {
                        list_no[k] = 0;
                    }
                    
           		}
            }
            //////////更新終わり/////////////
            //#pragma omp for
            for(i = 0; i < Number; i++)//力の計算
			{
                
                mlcl[i].f[X] = 0.0;
                mlcl[i].f[Y] = 0.0;
                mlcl[i].f[Z] = 0.0;
                
                /////////////////////////ボックスミュラー法によるノイズの生成/////////////////////////////
                p1 = sqrt(T*mlcl[i].r)*sqrt(-2.0*log(genrand_real3()));
				p2 = sqrt(T*mlcl[i].r)*sqrt(-2.0*log(genrand_real3()));
				theta = 2*PI*genrand_real2();
				psi = 2*PI*genrand_real2();
				mlcl[i].f[X] += p1*sin(theta)/sqrt(t);
				mlcl[i].f[Y] += p1*cos(theta)/sqrt(t);
				mlcl[i].f[Z] += p2*sin(psi)/sqrt(t);
                /////////////////////線形バネ接続//////////////////
                ////////ループの作成////////////
                //前側//
                mlcl[i].f[X] += LSx(i,mlcl[i].breakpair1_a,LINEAR);
                mlcl[i].f[Y] += LSy(i,mlcl[i].breakpair1_a,LINEAR);
                mlcl[i].f[Z] += LSz(i,mlcl[i].breakpair1_a,LINEAR);
                //後ろ側
                mlcl[i].f[X] += LSx(i,mlcl[i].breakpair1_b,LINEAR);
                mlcl[i].f[Y] += LSy(i,mlcl[i].breakpair1_b,LINEAR);
                mlcl[i].f[Z] += LSz(i,mlcl[i].breakpair1_b,LINEAR);
                
                ///////排除体積効果、衝突したらくっつく///////
                if(list_no[i] != 0)
                {
                    //printf("b");
                    for(j = 1; j < list_no[i] + 1; j++)
                    {
                        m = list_pair[i][j];
                        d = dist(i,m);
                        /////////くっつく//////////
                        if((mlcl[mlcl[i].bond_pair].bond_pair != i)&&(mlcl[i].bond_pair != -1))//片思いなら外す
                        {
                            mlcl[i].bond_pair = -1;
                        }
                        else if((d < mlcl[i].r + mlcl[m].r)//接触している
                                &&(m!=mlcl[i].breakpair1_a)&&(m!=mlcl[i].breakpair1_b)//前後の粒子でない
                                &&(mlcl[i].bond_pair == -1)&&(mlcl[m].bond_pair==-1)//くっつくペアを記憶していない
                                )
                        {
                            ////////くっつく先のペアを互いに記憶///////
                            mlcl[i].bond_pair = m;
                            mlcl[m].bond_pair = i;
                        }
                        
                        
                        
                        ///////くっついたもの同士の間に線形バネ//////////
                        if((mlcl[i].bond_pair == m)&&(mlcl[m].bond_pair == i))
                        {
                            
                            if(d>(mlcl[i].r+mlcl[m].r)*1.01)//確率で片思いに
                            {
                                mlcl[i].bond_pair = -1;
                            }
                            else//接続
                            {
                                mlcl[i].f[X] += LSx(i,m,LINEAR);
                                mlcl[i].f[Y] += LSy(i,m,LINEAR);
                                mlcl[i].f[Z] += LSz(i,m,LINEAR);
                                
                            }
                            
                        }
                        ///////////////////////////
                        ////////排除体積効果/////////
                        if((d<mlcl[i].r + mlcl[m].r)&&(m!=mlcl[i].breakpair1_a)&&(m!=mlcl[i].breakpair1_b))
                        {
                            mlcl[i].f[X] += EVEx(i, m, EXVE);
                            mlcl[i].f[Y] += EVEy(i, m, EXVE);
                            mlcl[i].f[Z] += EVEz(i, m, EXVE);
                        }
                        ///////////////////////////
                        if(d <= mlcl[i].r + mlcl[m].r)//衝突時にcount
                        {
                            //printf("b");
                            count[i][m] = count[i][m] + 1.0;
                        }
                    }
                }
                
                
                ///////////////////////////////////////////////////
                
				
                

			}
            for(i = 0; i < Number; i++)//新しい座標と速度に移行
            {
                //ノイズと排除体積効果と隣の粒子とのバネの力を付加
                //過減衰極限で計算
                mlcl[i].p_n[X] = mlcl[i].p[X] + t * mlcl[i].f[X]/mlcl[i].gamma;
				mlcl[i].p_n[Y] = mlcl[i].p[Y] + t * mlcl[i].f[Y]/mlcl[i].gamma;
				mlcl[i].p_n[Z] = mlcl[i].p[Z] + t * mlcl[i].f[Z]/mlcl[i].gamma;
                //位置の更新
                mlcl[i].p[X] = mlcl[i].p_n[X];
			    mlcl[i].p[Y] = mlcl[i].p_n[Y];
			    mlcl[i].p[Z] = mlcl[i].p_n[Z];

            }
            
            if(step%10000 == 0)
            {
                printf("\r%lld",step);
            }
            fflush(stdout);
            
            
        
            if((step%10000 == 0)&&(step!=0))
            {
                
                //テキストファイルへの書き出し
                FILE *fp;
                
                char filename[128];
                sprintf(filename, "resultfree/breakcountmatrix%lld.txt",step);
                fp = fopen(filename,"wt");
                
                //////////////////ファイル出力/////////////////////
                sum = 0.0;
                for(i = 0; i < Number; i++)
                {
                    for(j = 0; j < Number; j++)
                    {
                        sum += count[i][j];
                        if(Max<=count[i][j])
                        {
                            Max = count[i][j];
                        }
                        if(min>=count[i][j])
                        {
                            min = count[i][j];
                        }

                    }
                }
                for(i = 0; i < Number; i++)//衝突回数を行列表示
                {
                    for(j = 0; j < Number; j++)
                    {
                            fprintf(fp,"%lf ",count[i][j]/step);

                    }
                    fprintf(fp,"\n");

                }
                
                
                fclose(fp);
                /////////現状の保存//////////
                FILE *gp;
                char logfilename[128];
                sprintf(logfilename, "resultfree/position%lld.txt",step);
                gp = fopen(logfilename,"wt");
                fprintf(gp,"%d\n",large+small+large+small);
                fprintf(gp,"%d\n", large);
                fprintf(gp,"%d\n", small);
                fprintf(gp,"%lld\n", step);//step
                fprintf(gp,"%d\n", 2);//save
                for(i=0;i<Number;i++){
                    fprintf(gp, "%d %lf %lf %lf %lf %lf %lf\n", i, mlcl[i].p[X], mlcl[i].p[Y], mlcl[i].p[Z], mlcl[i].r, 1.0, mlcl[i].scale);
                }
                fclose(gp);
                ///////////////////////////
                
            }
            ///////////////出力おわり////////////////////
            
        }
        
	}
        


    //テキストファイルへの最終書き出し
    FILE *fp;
    fp = fopen("resultfree/breakcountmatrix.txt","wt");

    sum = 0;
    for(i = 0; i < Number; i++)
    {
        for(j = 0; j < Number; j++)
        {
            if((mlcl[i].scale == 1.0)&&(mlcl[j].scale == 1.0))
            {
                sum += count[i][j];
                if(Max<=count[i][j])
                {
                    Max = count[i][j];
                }
                if(min>=count[i][j])
                {
                    min = count[i][j];
                }
            }
        }
    }
    for(i = 0; i < Number; i++)//衝突回数を行列表示
    {
        for(j = 0; j < Number; j++)
        {
            fprintf(fp,"%lf ",count[i][j]/step);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}


int main(int argc,char** argv)
{
    time(&times);
    init_genrand(times);
	readposition(argc, argv);
    molsposition();
    free(mlcl);
	return(0);
}
