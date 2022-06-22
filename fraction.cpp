#include<iostream>
#include<fstream>
#include<cmath>
#define NP 70
#define NK 70
#define NMG 104
#define ZK 1
#define ZMG 2
using namespace std;

float strength[80][10000], I_avg[80], I_std[80];

struct example
{
    char name[5];
    float x, y, z;
};

example p_info[10000][80];
example k_info[10000][80];
//example mg_info[10000][200];

int main()
{
    ifstream infile;
    infile.open("p.pdb");

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

    infile.open("k.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将原子坐标格式化存入结构数组
    i = 1;
    while (!infile.eof())
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
        infile >> k_info[j][k].name >> k_info[j][k].x >> k_info[j][k].y >> k_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    /*infile.open("mg.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将原子坐标格式化存入结构数组
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
    infile.close();*/

    ofstream outfile, outfile1;
    outfile.open("fraction.dat", ios::out | ios::ate);
    outfile1.open("test.dat", ios::out | ios::ate);

    // 计算磷原子周围离子强度随时间的变化(K)
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= t; j++)
        {
            strength[i][j] = 0;
            for (k = 1; k <= NK; k++)
            {
                strength[i][j] += ZK / ((k_info[j][k].x - p_info[j][i].x) * (k_info[j][k].x - p_info[j][i].x) + (k_info[j][k].y - p_info[j][i].y) * (k_info[j][k].y - p_info[j][i].y) + (k_info[j][k].z - p_info[j][i].z) * (k_info[j][k].z - p_info[j][i].z));
            }
        }
    }

    for (i = 1; i <= t; i++)
    {
        outfile1 << i << " " << strength[7][i] << endl;
    }
    outfile1.close();

    /*// 计算磷原子周围离子强度随时间的变化(K+MG)
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= t; j++)
        {
            for (k = 1; k <= NMG; k++)
            {
                strength[i][j] += ZMG / ((mg_info[j][k].x - p_info[j][i].x) * (mg_info[j][k].x - p_info[j][i].x) + (mg_info[j][k].y - p_info[j][i].y) * (mg_info[j][k].y - p_info[j][i].y) + (mg_info[j][k].z - p_info[j][i].z) * (mg_info[j][k].z - p_info[j][i].z));
            }
        }
    }*/

    // 计算磷原子周围离子强度的平均值
    for (i = 1; i <= NP; i++)
    {
        I_avg[i] = 0;
        for (j = 1; j <= t; j++)
        {
            I_avg[i] += strength[i][j];
        }
        I_avg[i] = I_avg[i] / t;
    }

    // 计算磷原子周围离子强度的标准差
    for (i = 1; i <= NP; i++)
    {
        I_std[i] = 0;
        for (j = 1; j <= t; j++)
        {
            I_std[i] += (strength[i][j] - I_avg[i]) * (strength[i][j] - I_avg[i]);
        }
        I_std[i] = sqrt(I_std[i] / t);
    }

    // 打印结果
    for (i = 1; i <= NP; i++)
    {
        outfile << i << " " << I_avg[i] << " " << I_std[i] << " " << I_std[i] / I_avg[i] << endl;
    }
    outfile.close();

    return 0;
}
