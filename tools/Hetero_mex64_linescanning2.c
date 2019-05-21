#include <stdio.h>
#include <mex.h>
#include <math.h>
#include <stdlib.h>


//////////////////////////////////////////////////////////////
/// [Results,VLE_raw]=Hetero_mex(Volume,NumberOfPoints,Lmax);
//////////////////////// Input variable //////////////////////
// Volume : 3-Dimensional volume data (double type)
// NumberOfPoints : Number of pair of points
// Lmax : maximum number of length between two points
// Direction : generated unit vector data
////////////// Output variable ////////////
// Results : [L] [VLE_sum] [NumberOfPoints] [Mean]
// VLE_raw : Raw data
//////////////////////////////////////////////////////////////

double calc_mean(int* array, int size);
double calc_std(int* array, int size, double meanValue);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double* Volume=NULL;
    int *Dimension=NULL;
    double *NumberOfPoints=NULL;
    double *Direction=NULL;
    double unitVector[3]={0},HalfDist,distance;
    double *Results=NULL;
    double *VLE_raw=NULL;
    double c1,c2,c3;
    int NumberOfDirection=0;
    int NumOfLine,NumOfPlane=0;
    int Num_points;
    int LMax;
    int current_index;
    int p_dist;
    int i,j,k;
    int Dim_max,Dim_min;
    int cnt_VLE_raw=0;
    int VLE_raw_x,VLE_raw_y;
    
    Dimension=mxGetDimensions(prhs[0]);
    NumberOfPoints=mxGetPr(prhs[1]);
    Direction=mxGetPr(prhs[2]);
    NumberOfDirection=mxGetM(prhs[2]);
    Volume=mxGetPr( prhs[0] );
    
    c1=(Dimension[0]-1)*0.5;
    c2=(Dimension[1]-1)*0.5;
    c3=(Dimension[2]-1)*0.5;
    
    NumOfPlane=Dimension[0]*Dimension[1];
    NumOfLine=Dimension[0];
    
    
    /////////////////////////////////////////////////
    /* Calculate minimum and maximum of Dimensions */
    for (i=0;i<3;i++){
        Dim_max=Dimension[i];
        Dim_min=Dimension[i];
        if (Dimension[i] > Dim_max)
        {
            Dim_max = Dimension[i];
        }
        else if (Dimension[i] < Dim_min)
        {
            Dim_min = Dimension[i];
        }
    }
    
    LMax=(int)(Dim_min*0.5);
    /////////////////////////////////////////////////
    
    Num_points=(int)*NumberOfPoints;    
    plhs[0]=mxCreateDoubleMatrix(LMax,5,mxREAL);
    Results=mxGetPr(plhs[0]);
    
    /////////////////////////////////////////////////
    /* scanning with respect to the Length */
    for (i=0;i<LMax;i++){
//     for (i=5;i<=5;i++){
        double temp_mean1=0;
        double temp_mean2=0;
        int *VLE_vol_temp=NULL;      // 동적할당 방법 확인
        int *VLE_diff_temp=NULL;
        int num_vol_temp=0;
        int cnt_temp=1;
        
        num_vol_temp=(i+1)*(Num_points);
        
        VLE_diff_temp=(int(*))malloc(sizeof(int)*Num_points);
        VLE_vol_temp=(int(*))malloc(sizeof(int)*num_vol_temp);
        
        for (j=0;j<(int)*NumberOfPoints;j++){
            int p1,p2,p3,idx_direction;
            int PT11,PT12,PT13,PT21,PT22,PT23;
            int cnt_num_sum=0;
            double d1,d2,d3;
            double value1,value2,Dvalue;
            
            
            
            // rand() % (max_number + 1 - minimum_number)) + minimum_number
            p1=rand() % (Dimension[0] + 1 - 1) + 1;
            p2=rand() % (Dimension[1] + 1 - 1) + 1;
            p3=rand() % (Dimension[2] + 1 - 1) + 1;
            
            PT11=p1;    PT12=p2;    PT13=p3;
            
            idx_direction=rand() % (NumberOfDirection + 1 - 0) + 0;
            
            d1=Direction[idx_direction];
            d2=Direction[idx_direction+NumberOfDirection];
            d3=Direction[idx_direction+NumberOfDirection*2];
            
            PT21=p1+(int)(d1*i);    PT22=p2+(int)(d2*i);    PT23=p3+(int)(d3*i);
            
            // for the points out of boundary
            if ((PT11 >= Dimension[0])|(PT11 <= 0)|(PT12 >= Dimension[1])|(PT12 <= 0)|(PT13 >= Dimension[2])|(PT13 <= 0)|(PT21 >= Dimension[0])|(PT21 <= 0)|(PT22 >= Dimension[1])|(PT22 <= 0)|(PT23 >= Dimension[2])|(PT23 <= 0)){
                while ((PT11 >= Dimension[0])|(PT11 <= 0)|(PT12 >= Dimension[1])|(PT12 <= 0)|(PT13 >= Dimension[2])|(PT13 <= 0)|(PT21 >= Dimension[0])|(PT21 <= 0)|(PT22 >= Dimension[1])|(PT22 <= 0)|(PT23 >= Dimension[2])|(PT23 <= 0)){
                    p1=rand() % (Dimension[0] + 1 - 1) + 1;
                    p2=rand() % (Dimension[1] + 1 - 1) + 1;
                    p3=rand() % (Dimension[2] + 1 - 1) + 1;
                    PT11=p1;    PT12=p2;    PT13=p3;
                    PT21=p1+(int)(d1*i);    PT22=p2+(int)(d2*i);    PT23=p3+(int)(d3*i);
                }
            }
            
            for (k=0;k<=i;k++){
                int PT_crt1,PT_crt2,PT_crt3;
                int value_crt;
                
                PT_crt1=p1+(int)(d1*k);    PT_crt2=p2+(int)(d2*k);    PT_crt3=p3+(int)(d3*k);
                value_crt=Volume[PT_crt1+PT_crt2*NumOfLine+PT_crt3*NumOfPlane];
                cnt_temp++;
                VLE_vol_temp[k+j*(i+1)]=value_crt;
            }
            
            value1=Volume[PT11+PT12*NumOfLine+PT13*NumOfPlane];
            value2=Volume[PT21+PT22*NumOfLine+PT23*NumOfPlane];
            Dvalue=abs(value1-value2);
            VLE_diff_temp[j]=abs(value1-value2);
        }
        
        Results[i]=i;
        temp_mean1=calc_mean(VLE_diff_temp, Num_points);
        temp_mean2=calc_mean(VLE_vol_temp, num_vol_temp);
        Results[i+1*LMax]=temp_mean1;
        Results[i+2*LMax]=calc_std(VLE_diff_temp, Num_points,temp_mean1);
        Results[i+3*LMax]=temp_mean2;
        Results[i+4*LMax]=calc_std(VLE_vol_temp, num_vol_temp,temp_mean2);

        free(VLE_diff_temp);
        free(VLE_vol_temp);
    }
}



double calc_mean(int* array, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += array[i];
    }
    return sum / size;
}

double calc_std(int* array, int size, double meanValue , int option) {
    if (size < 2) return sqrt(-1.0);
    
    double sum = 0.0;
    double sd = 0.0;
    double diff;
    
    for (int i = 0; i < size; i++) {
        sum += (array[i] - meanValue) * (array[i] - meanValue);
    }
    sd = sqrt(sum / (size));
    return sd;
}
