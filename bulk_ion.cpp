#include<iostream>
#include<fstream>
#include<cmath>
#define NP 70
#define NK 70
#define NMG 104
#define dim 75
using namespace std;

int j_max[10000], k_max[10000];
float d_max_series[10000], d_temp, d_max, d_avg;
float center_jk_x[10000], center_jk_y[10000], center_jk_z[10000];
float N_bulk_K[10000], N_bulk_MG[10000], C_bulk_K[10000], C_bulk_MG[10000];

struct example
{
    char name[5];
    float x, y, z;
};

example p_info[10000][80];
example k_info[10000][80];
example mg_info[10000][110];

int main()
{
    ifstream file;
    file.open("p.pdb");

    if (!file)
    {
        cerr << "Open file failure!" << endl;
        return -1;
    }

    // 将P原子坐标格式化存入结构数组
    int i = 1, j, k, t;
    float aa;
    while (!file.eof())
    {
        aa = (i + 0.0) / (NP + 0.0); 
        k = fmod(i, NP);
        if (k == 0) 
        {
            k = NP;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        file >> p_info[j][k].name >> p_info[j][k].x >> p_info[j][k].y >> p_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    file.close();

    file.open("k.pdb");

    if (!file)
    {
        cerr << "Open file failure!" << endl;
        return -1;
    }

    // 将K离子坐标格式化存入结构数组
    i = 1;
    while (!file.eof())
    {
        aa = (i + 0.0) / (NK + 0.0); 
        k = fmod(i, NK);
        if (k == 0) 
        {
            k = NK;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        file >> k_info[j][k].name >> k_info[j][k].x >> k_info[j][k].y >> k_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    file.close();

    file.open("mg.pdb");

    if (!file)
    {
        cerr << "Open file failure!" << endl;
        return -1;
    }

    // 将MG离子坐标格式化存入结构数组
    i = 1;
    while (!file.eof())
    {
        aa = (i + 0.0) / (NMG + 0.0); 
        k = fmod(i, NMG);
        if (k == 0) 
        {
            k = NMG;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        file >> mg_info[j][k].name >> mg_info[j][k].x >> mg_info[j][k].y >> mg_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    file.close();

    ofstream outfile;
    outfile.open("bulk_ion.dat", ios::out | ios::ate);

    // 计算核酸分子中两两P原子之间最大距离随时间的变化，并记录P原子序号
    for (i = 1; i <= t; i++)
    {
        d_max_series[i] = 0;
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= NP; k++)
            {
                d_temp = sqrt((p_info[i][j].x - p_info[i][k].x) * (p_info[i][j].x - p_info[i][k].x) + (p_info[i][j].y - p_info[i][k].y) * (p_info[i][j].y - p_info[i][k].y) + (p_info[i][j].z - p_info[i][k].z) * (p_info[i][j].z - p_info[i][k].z));
                if (d_temp > d_max_series[i])
                {
                    d_max_series[i] = d_temp;
                    j_max[i] = j;
                    k_max[i] = k;
                }
            }
        }
    }

    // 计算两两P原子中点坐标随时间的变化
    for (i = 1; i <= t; i++)
    {
        center_jk_x[i] = (p_info[i][j_max[i]].x + p_info[i][k_max[i]].x) / 2;
        center_jk_y[i] = (p_info[i][j_max[i]].y + p_info[i][k_max[i]].y) / 2;
        center_jk_z[i] = (p_info[i][j_max[i]].z + p_info[i][k_max[i]].z) / 2;
    }

    // 计算核酸分子中两两P原子之间最大距离随时间变化的最大值
    d_max = 0;
    for (i = 1; i <= t; i++)
    {
        if (d_max_series[i] > d_max)
            d_max = d_max_series[i];
    }
    cout << d_max << endl;

    // 计算核酸分子中两两P原子之间最大距离随时间变化的平均值
    d_avg = 0;
    for (i = 1; i <= t; i++)
    {
        d_avg += d_max_series[i];
    }
    d_avg /= t;
    cout << d_avg << endl;

    // 计算K离子体浓度随时间的变化
    for (i = 1; i <= t; i++)
    {
        N_bulk_K[i] = 0;
        for (j = 1; j <= NK; j++)
        {
            d_temp = sqrt((k_info[i][j].x - center_jk_x[i]) * (k_info[i][j].x - center_jk_x[i]) + (k_info[i][j].y - center_jk_y[i]) * (k_info[i][j].y - center_jk_y[i]) + (k_info[i][j].z - center_jk_z[i]) * (k_info[i][j].z - center_jk_z[i]));
            if (d_temp > dim / 2 && d_temp < dim / 2 + 10)
                N_bulk_K[i]++;
        }
        C_bulk_K[i] = 1000 * N_bulk_K[i] / 0.000602 / 4.19 / (pow(dim / 2 + 10, 3) - pow(dim / 2, 3));
    }

    // 计算MG离子体浓度随时间的变化
    for (i = 1; i <= t; i++)
    {
        N_bulk_MG[i] = 0;
        for (j = 1; j <= NMG; j++)
        {
            d_temp = sqrt((mg_info[i][j].x - center_jk_x[i]) * (mg_info[i][j].x - center_jk_x[i]) + (mg_info[i][j].y - center_jk_y[i]) * (mg_info[i][j].y - center_jk_y[i]) + (mg_info[i][j].z - center_jk_z[i]) * (mg_info[i][j].z - center_jk_z[i]));
            if (d_temp > dim / 2 && d_temp < dim / 2 + 10)
                N_bulk_MG[i]++;
        }
        C_bulk_MG[i] = 1000 * N_bulk_MG[i] / 0.000602 / 4.19 / (pow(dim / 2 + 10, 3) - pow(dim / 2, 3));
        outfile << i << " " << d_max_series[i] << " " << C_bulk_K[i] << " " << C_bulk_MG[i] << endl;
    }
    outfile.close();

    return 0;
}
