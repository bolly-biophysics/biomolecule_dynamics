#include<iostream>
#include<fstream>
#include<cmath>
#define N0 2257
using namespace std;

float x_10_center[10000], y_10_center[10000], z_10_center[10000], x_35_center[10000], y_35_center[10000], z_35_center[10000];
float v1_10_x[10000], v1_10_y[10000], v1_10_z[10000], v2_10_x[10000], v2_10_y[10000], v2_10_z[10000];
float v1_35_x[10000], v1_35_y[10000], v1_35_z[10000], v2_35_x[10000], v2_35_y[10000], v2_35_z[10000];
float n_10_x[10000], n_10_y[10000], n_10_z[10000], n_35_x[10000], n_35_y[10000], n_35_z[10000];
float r_10_35_x[10000], r_10_35_y[10000], r_10_35_z[10000];
float norm_n_10[10000], norm_n_35[10000], norm_r_10_35[10000];
float theta[10000], tau[10000], d[10000];

struct example_info
{
    char name[5];
    float x, y, z;
};

example_info nucleic_info[10000][2300];

int main()
{
    ifstream infile;
    infile.open("nucleic.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将原子坐标格式化存入结构数组
    int i = 1, j, k, t;
    float aa;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (N0 + 0.0); 
        k = fmod(i, N0);
        if (k == 0) 
        {
            k = N0;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> nucleic_info[j][k].name >> nucleic_info[j][k].x >> nucleic_info[j][k].y >> nucleic_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    ofstream outfile;
    outfile.open("stacking.dat", ios::out | ios::ate);

    // 求碱基U10与U35芳香环的中心坐标(即N3与C6的中点)
    for (i = 1; i <= t; i++)
    {
        x_10_center[i] = (nucleic_info[i][295].x + nucleic_info[i][301].x) / 2;
        y_10_center[i] = (nucleic_info[i][295].y + nucleic_info[i][301].y) / 2;
        z_10_center[i] = (nucleic_info[i][295].z + nucleic_info[i][301].z) / 2;

        x_35_center[i] = (nucleic_info[i][1096].x + nucleic_info[i][1102].x) / 2;
        y_35_center[i] = (nucleic_info[i][1096].y + nucleic_info[i][1102].y) / 2;
        z_35_center[i] = (nucleic_info[i][1096].z + nucleic_info[i][1102].z) / 2;
    }

    // 分别求出U10与U35芳香环平面上的两个向量(C6-center与C5-center)
    for (i = 1; i <= t; i++)
    {
        v1_10_x[i] = nucleic_info[i][295].x - x_10_center[i];
        v1_10_y[i] = nucleic_info[i][295].y - y_10_center[i];
        v1_10_z[i] = nucleic_info[i][295].z - z_10_center[i];

        v2_10_x[i] = nucleic_info[i][297].x - x_10_center[i];
        v2_10_y[i] = nucleic_info[i][297].y - y_10_center[i];
        v2_10_z[i] = nucleic_info[i][297].z - z_10_center[i];

        v1_35_x[i] = nucleic_info[i][1096].x - x_35_center[i];
        v1_35_y[i] = nucleic_info[i][1096].y - y_35_center[i];
        v1_35_z[i] = nucleic_info[i][1096].z - z_35_center[i];

        v2_35_x[i] = nucleic_info[i][1098].x - x_35_center[i];
        v2_35_y[i] = nucleic_info[i][1098].y - y_35_center[i];
        v2_35_z[i] = nucleic_info[i][1098].z - z_35_center[i];
    }

    // 分别求出U10与U35芳香环平面的法向量n_10与n_35
    for (i = 1; i <= t; i++)
    {
        n_10_x[i] = v1_10_y[i] * v2_10_z[i] - v1_10_z[i] * v2_10_y[i];
        n_10_y[i] = v1_10_z[i] * v2_10_x[i] - v1_10_x[i] * v2_10_z[i];
        n_10_z[i] = v1_10_x[i] * v2_10_y[i] - v1_10_y[i] * v2_10_x[i];
        norm_n_10[i] = sqrt(n_10_x[i] * n_10_x[i] + n_10_y[i] * n_10_y[i] + n_10_z[i] * n_10_z[i]);

        n_35_x[i] = v1_35_y[i] * v2_35_z[i] - v1_35_z[i] * v2_35_y[i];
        n_35_y[i] = v1_35_z[i] * v2_35_x[i] - v1_35_x[i] * v2_35_z[i];
        n_35_z[i] = v1_35_x[i] * v2_35_y[i] - v1_35_y[i] * v2_35_x[i];
        norm_n_35[i] = sqrt(n_35_x[i] * n_35_x[i] + n_35_y[i] * n_35_y[i] + n_35_z[i] * n_35_z[i]);
    }

    // 计算n_10与n_35之间的夹角
    for (i = 1; i <= t; i++)
    {
        theta[i] = acos((n_10_x[i] * n_35_x[i] + n_10_y[i] * n_35_y[i] + n_10_z[i] * n_35_z[i]) / norm_n_10[i] / norm_n_35[i]) / 3.14 * 180;
    }

    // 求U10与U35芳香环平面中心连线的向量
    for (i = 1; i <= t; i++)
    {
        r_10_35_x[i] = x_10_center[i] - x_35_center[i];
        r_10_35_y[i] = y_10_center[i] - y_35_center[i];
        r_10_35_z[i] = z_10_center[i] - z_35_center[i];
        norm_r_10_35[i] = sqrt(r_10_35_x[i] * r_10_35_x[i] + r_10_35_y[i] * r_10_35_y[i] + r_10_35_z[i] * r_10_35_z[i]);
    }

    // 计算r_10_35与n_35之间的夹角
    for (i = 1; i <= t; i++)
    {
        tau[i] = acos((r_10_35_x[i] * n_35_x[i] + r_10_35_y[i] * n_35_y[i] + r_10_35_z[i] * n_35_z[i]) / norm_r_10_35[i] / norm_n_35[i]) / 3.14 * 180;
    }

    // 计算r_10_35在n_35上的投影长度
    for (i = 1; i <= t; i++)
    {
        d[i] = norm_r_10_35[i] * cos(tau[i] / 180 * 3.14);
        outfile << i << " " << theta[i] << " " << tau[i] << " " << norm_r_10_35[i] << endl;
    }
    outfile << endl;

    return 0;
}
