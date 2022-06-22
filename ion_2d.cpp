#include<iostream>
#include<fstream>
#include<cmath>
#define NP 70
#define NK 70
#define NMG 104
#define R 8
using namespace std;

float density_k[80][80][10000], density_k_avg[80][80], density_k_std[80][80];
float density_mg[80][80][10000], density_mg_avg[80][80], density_mg_std[80][80];
float I_k[80][80][10000], I_k_avg[80][80], I_k_std[80][80];
float I_mg[80][80][10000], I_mg_avg[80][80], I_mg_std[80][80];

struct example
{
    char name[5];
    float x, y, z;
};

example p_info[10000][80];
example k_info[10000][80];
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
    int i = 1, j, k, l, t;
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

    // 将K离子坐标格式化存入结构数组
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

    infile.open("mg.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将MG离子坐标格式化存入结构数组
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

    ofstream outfile, outfile1, outfile2, outfile3;
    outfile.open("k_2d.dat", ios::out | ios::ate);
    outfile1.open("mg_2d.dat", ios::out | ios::ate);
    outfile2.open("I_k_2d.dat", ios::out | ios::ate);
    outfile3.open("I_mg_2d.dat", ios::out | ios::ate);

    // 计算两两磷原子周围K离子个数随时间的变化
    float d1, d2;
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                density_k[i][j][k] = 0;
                for (l = 1; l <= NK; l++)
                {
                    d1 = sqrt((k_info[k][l].x - p_info[k][i].x) * (k_info[k][l].x - p_info[k][i].x) + (k_info[k][l].y - p_info[k][i].y) * (k_info[k][l].y - p_info[k][i].y) + (k_info[k][l].z - p_info[k][i].z) * (k_info[k][l].z - p_info[k][i].z));
                    d2 = sqrt((k_info[k][l].x - p_info[k][j].x) * (k_info[k][l].x - p_info[k][j].x) + (k_info[k][l].y - p_info[k][j].y) * (k_info[k][l].y - p_info[k][j].y) + (k_info[k][l].z - p_info[k][j].z) * (k_info[k][l].z - p_info[k][j].z));
                    if (d1 < R && d2 < R)
                        density_k[i][j][k]++;
                }
            }
        }
    }

    // 计算两两磷原子周围K离子个数的平均值
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            density_k_avg[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_k_avg[i][j] += density_k[i][j][k];
            }
            density_k_avg[i][j] = density_k_avg[i][j] / t;
        }
    }

    // 计算两两磷原子周围K离子个数的标准差
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {        
            density_k_std[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_k_std[i][j] += (density_k[i][j][k] - density_k_avg[i][j]) * (density_k[i][j][k] - density_k_avg[i][j]);
            }
            density_k_std[i][j] = sqrt(density_k_std[i][j] / t);
        }
    }

    // 计算两两磷原子周围MG离子个数随时间的变化
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                density_mg[i][j][k] = 0;
                for (l = 1; l <= NMG; l++)
                {
                    d1 = sqrt((mg_info[k][l].x - p_info[k][i].x) * (mg_info[k][l].x - p_info[k][i].x) + (mg_info[k][l].y - p_info[k][i].y) * (mg_info[k][l].y - p_info[k][i].y) + (mg_info[k][l].z - p_info[k][i].z) * (mg_info[k][l].z - p_info[k][i].z));
                    d2 = sqrt((mg_info[k][l].x - p_info[k][j].x) * (mg_info[k][l].x - p_info[k][j].x) + (mg_info[k][l].y - p_info[k][j].y) * (mg_info[k][l].y - p_info[k][j].y) + (mg_info[k][l].z - p_info[k][j].z) * (mg_info[k][l].z - p_info[k][j].z));
                    if (d1 < R && d2 < R)
                        density_mg[i][j][k]++;
                }
            }
        }
    }

    // 计算两两磷原子周围MG离子个数的平均值
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            density_mg_avg[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_mg_avg[i][j] += density_mg[i][j][k];
            }
            density_mg_avg[i][j] = density_mg_avg[i][j] / t;
        }
    }

    // 计算两两磷原子周围MG离子个数的标准差
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {        
            density_mg_std[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_mg_std[i][j] += (density_mg[i][j][k] - density_mg_avg[i][j]) * (density_mg[i][j][k] - density_mg_avg[i][j]);
            }
            density_mg_std[i][j] = sqrt(density_mg_std[i][j] / t);
        }
    }

    // 输出结果
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            outfile << i << " " << j << " " << density_k_avg[i][j] << " " << density_k_std[i][j] << endl;
            outfile1 << i << " " << j << " " << density_mg_avg[i][j] << " " << density_mg_std[i][j] << endl;
        }
        outfile << endl;
        outfile1 << endl;
    }
    outfile.close();
    outfile1.close();

    // 计算两两磷原子周围K离子强度随时间的变化
    float square_d1, square_d2;
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                I_k[i][j][k] = 0;
                for (l = 1; l <= NK; l++)
                {
                    square_d1 = (k_info[k][l].x - p_info[k][i].x) * (k_info[k][l].x - p_info[k][i].x) + (k_info[k][l].y - p_info[k][i].y) * (k_info[k][l].y - p_info[k][i].y) + (k_info[k][l].z - p_info[k][i].z) * (k_info[k][l].z - p_info[k][i].z);
                    square_d2 = (k_info[k][l].x - p_info[k][j].x) * (k_info[k][l].x - p_info[k][j].x) + (k_info[k][l].y - p_info[k][j].y) * (k_info[k][l].y - p_info[k][j].y) + (k_info[k][l].z - p_info[k][j].z) * (k_info[k][l].z - p_info[k][j].z);
                    I_k[i][j][k] += 1 / (square_d1 + square_d2);
                }
            }
        }
    }

    // 计算两两磷原子周围K离子强度的平均值
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            I_k_avg[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                I_k_avg[i][j] += I_k[i][j][k];
            }
            I_k_avg[i][j] = I_k_avg[i][j] / t;
        }
    }

    // 计算两两磷原子周围K离子强度的标准差
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {        
            I_k_std[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                I_k_std[i][j] += (I_k[i][j][k] - I_k_avg[i][j]) * (I_k[i][j][k] - I_k_avg[i][j]);
            }
            I_k_std[i][j] = sqrt(I_k_std[i][j] / t);
        }
    }

    // 计算两两磷原子周围MG离子强度随时间的变化
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                I_mg[i][j][k] = 0;
                for (l = 1; l <= NMG; l++)
                {
                    square_d1 = (mg_info[k][l].x - p_info[k][i].x) * (mg_info[k][l].x - p_info[k][i].x) + (mg_info[k][l].y - p_info[k][i].y) * (mg_info[k][l].y - p_info[k][i].y) + (mg_info[k][l].z - p_info[k][i].z) * (mg_info[k][l].z - p_info[k][i].z);
                    square_d2 = (mg_info[k][l].x - p_info[k][j].x) * (mg_info[k][l].x - p_info[k][j].x) + (mg_info[k][l].y - p_info[k][j].y) * (mg_info[k][l].y - p_info[k][j].y) + (mg_info[k][l].z - p_info[k][j].z) * (mg_info[k][l].z - p_info[k][j].z);
                    I_mg[i][j][k] += 2 / (square_d1 + square_d2);
                }
            }
        }
    }

    // 计算两两磷原子周围MG离子强度的平均值
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            I_mg_avg[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                I_mg_avg[i][j] += I_mg[i][j][k];
            }
            I_mg_avg[i][j] = I_mg_avg[i][j] / t;
        }
    }

    // 计算两两磷原子周围MG离子强度的标准差
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {        
            I_mg_std[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                I_mg_std[i][j] += (I_mg[i][j][k] - I_mg_avg[i][j]) * (I_mg[i][j][k] - I_mg_avg[i][j]);
            }
            I_mg_std[i][j] = sqrt(I_mg_std[i][j] / t);
        }
    }

    // 输出结果
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            outfile2 << i << " " << j << " " << I_k_avg[i][j] << " " << I_k_std[i][j] << endl;
            outfile3 << i << " " << j << " " << I_mg_avg[i][j] << " " << I_mg_std[i][j] << endl;
        }
        outfile2 << endl;
        outfile3 << endl;
    }
    outfile2.close();
    outfile3.close();

    return 0;
}
