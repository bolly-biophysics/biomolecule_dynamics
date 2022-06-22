#include<iostream>
#include<fstream>
#include<cmath>
#define NP 70
#define NMG 104
using namespace std;

float fraction_I[80][10000], fraction_I_avg[80], fraction_I_std[80];

struct example
{
    char name[5];
    float x, y, z;
};

example p_info[10000][80];
example mg_info[10000][200];

int main()
{
    ifstream infile;
    infile.open("p.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将P原子坐标格式化存入结构数组
    int i = 1, j, k, t;
    float aa;
    while (!infile.eof())
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
        infile >> p_info[j][k].name >> p_info[j][k].x >> p_info[j][k].y >> p_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    infile.open("mg.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将K离子坐标格式化存入结构数组
    i = 1;
    while (!infile.eof())
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
        infile >> mg_info[j][k].name >> mg_info[j][k].x >> mg_info[j][k].y >> mg_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    ofstream outfile, outfile1;
    outfile.open("mg_ep_series.dat", ios::out | ios::ate);
    outfile1.open("mg_ep.dat", ios::out | ios::ate);

    // 计算每个磷原子周围K离子强度随时间的变化
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= t; j++)
        {
            fraction_I[i][j] = 0;
            for (k = 1; k <= NMG; k++)
            {
                fraction_I[i][j] += 2 / sqrt((mg_info[j][k].x - p_info[j][i].x) * (mg_info[j][k].x - p_info[j][i].x) + (mg_info[j][k].y - p_info[j][i].y) * (mg_info[j][k].y - p_info[j][i].y) + (mg_info[j][k].z - p_info[j][i].z) * (mg_info[j][k].z - p_info[j][i].z));
            }
        }
    }

    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            outfile << i << " " << j << " " << fraction_I[j][i] << endl;
        }
        outfile << endl;
    }
    outfile.close();

    // 计算每个磷原子周围K离子密度的平均值
    for (i = 1; i <= NP; i++)
    {
        fraction_I_avg[i] = 0;
        for (j = 1; j <= t; j++)
        {
            fraction_I_avg[i] += fraction_I[i][j];
        }
        fraction_I_avg[i] /= t;
    }

    // 计算每个磷原子周围K离子密度的标准差
    for (i = 1; i <= NP; i++)
    {
        fraction_I_std[i] = 0;
        for (j = 1; j <= t; j++)
        {        
            fraction_I_std[i] += (fraction_I[i][j] - fraction_I_avg[i]) * (fraction_I[i][j] - fraction_I_avg[i]);
        }
        fraction_I_std[i] = sqrt(fraction_I_std[i] / t);
        outfile1 << i << " " << fraction_I_avg[i] << " " << fraction_I_std[i] / fraction_I_avg[i] << endl;
    }
    outfile1.close();

    return 0;
}
