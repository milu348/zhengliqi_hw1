#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if(a.rows!=b.rows||a.cols!=b.cols){
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }else{
        Matrix c = create_matrix(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
            {
                c.data[i][j] = a.data[i][j] + b.data[i][j];
            }
        }
        return c;

    }
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if(a.rows!=b.rows||a.cols!=b.cols){
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }else{
        Matrix c = create_matrix(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
            {
                c.data[i][j] = a.data[i][j] - b.data[i][j];
            }
        }
        return c;

    }
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if(a.cols != b.rows){
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0,0);
    }else{
        Matrix c = create_matrix(a.rows,b.cols);
        for(int i = 0;i < a.rows;i++){
            for(int j = 0;j < b.cols; j++){
                c.data[i][j] = 0;                 //初始化
                for(int k =0;k < a.cols;k++){
                    c.data[i][j] += a.data[i][k]*b.data[k][j];
                }
            }
        } 
        return c;
    }
}

Matrix scale_matrix(Matrix a, double k)
{
    for(int i = 0;i < a.rows;i++){
        for(int j = 0;j < a.cols;j++){
            a.data[i][j] *= k;
        }
    }
    return a;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix b = create_matrix(a.cols,a.rows);
    for(int i = 0;i < a.cols;i++){
        for(int j = 0;j < a.rows;j++){
            b.data[i][j] = a.data[j][i];
        }
    }
    return b;
}

double det_matrix(Matrix a)
{
    if(a.rows != a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }else{                                                  //拉普拉斯定理+递归
        if(a.rows == 0){
            return 0;
        }else if(a.rows == 1){
            return a.data[0][0];
        }else{
            double d = 0;
            for(int i = 0;i < a.rows;i++){
                Matrix b = create_matrix(a.rows - 1,a.cols - 1); //余子式
                for(int j = 0;j < b.rows;j++){
                    for(int k = 0;k < b.cols;k++){
                        if(k < i){
                            b.data[j][k] = a.data[j+1][k];
                        }else{
                            b.data[j][k] = a.data[j+1][k+1];
                        }
                    }
                }
                if(b.rows == 1){
                    d += a.data[0][i]*b.data[0][0]*pow(-1,i); 
                }else{
                d += a.data[0][i]*det_matrix(b)*pow(-1,i); 
                }
            }
            return d; 
        }
    }
}

Matrix inv_matrix(Matrix a)
{
    if(a.rows != a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0,0);
    }else if(det_matrix(a) == 0){
        printf("Error: The matrix is singular.\n");
        return create_matrix(0,0);
    }else{
        double det_a = det_matrix(a);
        Matrix b = create_matrix(a.rows,a.cols);
        for(int i = 0;i < a.rows;i++){
            for(int j = 0;j < a.cols;j++){
                Matrix c = create_matrix(a.rows - 1,a.cols - 1);
                int flag_k = 0;
                int flag_l = 0;
                for(int k = 0;k < c.rows;k++){            //a.data[i][j]的余子式
                    for(int l = 0;l < c.cols;l++){
                        if(k == i){
                            flag_k = 1;
                        }
                        if(l == j){
                            flag_l = 1;
                        }else{
                            flag_l = 0;                    //k不会出现与i比较变化，但l在不同k时会重置，要重新和j比较
                        }
                        c.data[k][l] = a.data[k+flag_k][l+flag_l];
                    }
                }
                b.data[j][i] = det_matrix(c)*pow(-1,i+j)/det_a;
            }
        }
        return b;
    }
}

int rank_matrix(Matrix a)
{
    int min = a.rows < a.cols ? a.rows : a.cols; //行列数最小值
    int rank = 0;
    int flag = 0;
    int col = 0;
    for(int i = 0;i < min;i++){
        if(a.data[i-col][i] != 0){                   //对角线元素不为0
            rank++;
            for(int j = i-col+1;j < a.rows;j++){
                if(a.data[j][i] != 0){           
                    double k = a.data[j][i]/a.data[i][i]; //消元系数
                    for(int l = i;l < a.cols;l++){
                        a.data[j][l] -= k*a.data[i-col][l];
                    }
                }
            }
        }else{
            flag = 0;
            for(int j = i-col+1;j < a.rows;j++){
                if(a.data[j][i] != 0){               //对角线元素为0，寻找非零行
                    for(int k = i;k < a.cols;k++){
                        double temp = a.data[i-col][k];
                        a.data[i-col][k] = a.data[j][k];
                        a.data[j][k] = temp;
                    }
                    flag=1;
                    break;
                }
            }
            if(flag == 1){
                rank++;
                for(int j = i-col+1;j < a.rows;j++){
                    if(a.data[j][i] != 0){           
                        double k = a.data[j][i]/a.data[i][i]; //消元系数
                        for(int l = i;l < a.cols;l++){
                            a.data[j][l] -= k*a.data[i-col][l];
                        }
                    }
                }
            }else{
                col++;
            }
        }
    }
    return rank;
}

double trace_matrix(Matrix a)
{
    if(a.rows != a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }else{
        double t = 0;
        for(int i = 0;i < a.rows;i++){
            t += a.data[i][i];
        }
        return t;
    }
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}