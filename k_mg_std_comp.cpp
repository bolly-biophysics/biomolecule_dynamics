#include<iostream>
#include<fstream>
#include<cmath>
#define NP 70
using namespace std;

float d_k[80][80][10000], d_mg[80][80][10000], avg_d_k[80][80], std_d_k[80][80], avg_d_mg[80][80], std_d_mg[80][80];
float ang_k[80][80][10000], ang_mg[80][80][10000], avg_ang_k[80][80], std_ang_k[80][80], avg_ang_mg[80][80], std_ang_mg[80][80];
float vector_k_x[80][10000], vector_k_y[80][10000], vector_k_z[80][10000];
float vector_mg_x[80][10000], vector_mg_y[80][10000], vector_mg_z[80][10000];

struct example_info
{
    char name[5];
    float x, y, z;
};

example_info p_info_k[10000][80];
example_info p_info_mg[10000][80];

int main()
{
    ifstream infile1, infile2;
    infile1.open("p_k.pdb");
    infile2.open("p_mg.pdb");

    if (!infile1)
    {
        cerr << "Open infile1 failure!" << endl;
        return -1;
    }

    // 将磷原子坐标格式化存入结构数组(K+)
    int i = 1, j, k, t;
    float aa;
    while (!infile1.eof())
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
        infile1 >> p_info_k[j][k].name >> p_info_k[j][k].x >> p_info_k[j][k].y >> p_info_k[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile1.close();

    if (!infile2)
    {
        cerr << "Open infile2 failure!" << endl;
        return -1;
    }

    // 将磷原子坐标格式化存入结构数组(MG)
    i = 1;
    while (!infile2.eof())
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
        infile2 >> p_info_mg[j][k].name >> p_info_mg[j][k].x >> p_info_mg[j][k].y >> p_info_mg[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile2.close();

    ofstream outfile1, outfile2, outfile3, outfile4;
    outfile1.open("std_d_k.dat", ios::out | ios::ate);
    outfile2.open("std_d_mg.dat", ios::out | ios::ate);
    outfile3.open("std_ang_k.dat", ios::out | ios::ate);
    outfile4.open("std_ang_mg.dat", ios::out | ios::ate);

    // 计算两两磷原子之间距离随时间的变化(K+)
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                d_k[i][j][k] = sqrt((p_info_k[k][i].x - p_info_k[k][j].x) * (p_info_k[k][i].x - p_info_k[k][j].x) + (p_info_k[k][i].y - p_info_k[k][j].y) * (p_info_k[k][i].y - p_info_k[k][j].y) + (p_info_k[k][i].z - p_info_k[k][j].z) * (p_info_k[k][i].z - p_info_k[k][j].z));
            }
        }
    }

    // 计算两两磷原子之间距离的平均值(K+)
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                avg_d_k[i][j] += d_k[i][j][k];
            }
            avg_d_k[i][j] = avg_d_k[i][j] / t;
        }
    }

    // 计算两两磷原子之间距离的标准差(K+)
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                std_d_k[i][j] += (d_k[i][j][k] - avg_d_k[i][j]) * (d_k[i][j][k] - avg_d_k[i][j]);
            }
            std_d_k[i][j] = sqrt(std_d_k[i][j] / t);
        }
    }

    // 计算两两磷原子之间距离随时间的变化(MG)
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                d_mg[i][j][k] = sqrt((p_info_mg[k][i].x - p_info_mg[k][j].x) * (p_info_mg[k][i].x - p_info_mg[k][j].x) + (p_info_mg[k][i].y - p_info_mg[k][j].y) * (p_info_mg[k][i].y - p_info_mg[k][j].y) + (p_info_mg[k][i].z - p_info_mg[k][j].z) * (p_info_mg[k][i].z - p_info_mg[k][j].z));
            }
        }
    }

    // 计算两两磷原子之间距离的平均值(MG)
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                avg_d_mg[i][j] += d_mg[i][j][k];
            }
            avg_d_mg[i][j] = avg_d_mg[i][j] / t;
        }
    }

    // 计算两两磷原子之间距离的标准差(MG)
    for (i = 1; i <= NP; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            for (k = 1; k <= t; k++)
            {
                std_d_mg[i][j] += (d_mg[i][j][k] - avg_d_mg[i][j]) * (d_mg[i][j][k] - avg_d_mg[i][j]);
            }
            std_d_mg[i][j] = sqrt(std_d_mg[i][j] / t);
        }
    }

    // 打印结果
    for (i = 1; i <= NP; i++)
        {
            for (j = 1; j <= NP; j++)
            {        
                outfile1 << i << " " << j << " " << std_d_k[i][j] << endl;
                outfile2 << i << " " << j << " " << std_d_mg[i][j] << endl;
            }
        outfile1 << endl;
        outfile2 << endl;
        }
    outfile1.close();
    outfile2.close();

    // 计算相邻磷原子之间形成的矢量随时间的变化(K+)
    for (i = 1; i <= NP - 1; i++)
    {
        for (j = 1; j <= t; j++)
        {
            vector_k_x[i][j] = p_info_k[j][i + 1].x - p_info_k[j][i].x;
            vector_k_y[i][j] = p_info_k[j][i + 1].y - p_info_k[j][i].y;
            vector_k_z[i][j] = p_info_k[j][i + 1].z - p_info_k[j][i].z;
        }
    }

    // 计算两两矢量之间夹角余弦随时间的变化(K+)
    for (i = 1; i <= NP - 1; i++)
    {
        for (j = 1; j <= NP - 1; j++)
        {
            for (k = 1; k <= t; k++)
            {
                ang_k[i][j][k] = (vector_k_x[i][k] * vector_k_x[j][k] + vector_k_y[i][k] * vector_k_y[j][k] + vector_k_z[i][k] * vector_k_z[j][k]) / sqrt(vector_k_x[i][k] * vector_k_x[i][k] + vector_k_y[i][k] * vector_k_y[i][k] + vector_k_z[i][k] * vector_k_z[i][k]) / sqrt(vector_k_x[j][k] * vector_k_x[j][k] + vector_k_y[j][k] * vector_k_y[j][k] + vector_k_z[j][k] * vector_k_z[j][k]);
            }
        }
    }

    // 计算余弦的平均值(K+)
    for (i = 1; i <= NP - 1; i++)
    {
        for (j = 1; j <= NP - 1; j++)
        {
            for (k = 1; k <= t; k++)
            {
                avg_ang_k[i][j] += ang_k[i][j][k];
            }
            avg_ang_k[i][j] = avg_ang_k[i][j] / t;
        }
    }

    // 计算余弦的标准差(K+)
    for (i = 1; i <= NP - 1; i++)
    {
        for (j = 1; j <= NP - 1; j++)
        {
            for (k = 1; k <= t; k++)
            {
                std_ang_k[i][j] += (ang_k[i][j][k] - avg_ang_k[i][j]) * (ang_k[i][j][k] - avg_ang_k[i][j]);
            }
            std_ang_k[i][j] = sqrt(std_ang_k[i][j] / t);
        }
    }

    // 计算相邻磷原子之间形成的矢量随时间的变化(MG)
    for (i = 1; i <= NP - 1; i++)
    {
        for (j = 1; j <= t; j++)
        {
            vector_mg_x[i][j] = p_info_mg[j][i + 1].x - p_info_mg[j][i].x;
            vector_mg_y[i][j] = p_info_mg[j][i + 1].y - p_info_mg[j][i].y;
            vector_mg_z[i][j] = p_info_mg[j][i + 1].z - p_info_mg[j][i].z;
        }
    }

    // 计算两两矢量之间夹角余弦随时间的变化(MG)
    for (i = 1; i <= NP - 1; i++)
    {
        for (j = 1; j <= NP - 1; j++)
        {
            for (k = 1; k <= t; k++)
            {
                ang_mg[i][j][k] = (vector_mg_x[i][k] * vector_mg_x[j][k] + vector_mg_y[i][k] * vector_mg_y[j][k] + vector_mg_z[i][k] * vector_mg_z[j][k]) / sqrt(vector_mg_x[i][k] * vector_mg_x[i][k] + vector_mg_y[i][k] * vector_mg_y[i][k] + vector_mg_z[i][k] * vector_mg_z[i][k]) / sqrt(vector_mg_x[j][k] * vector_mg_x[j][k] + vector_mg_y[j][k] * vector_mg_y[j][k] + vector_mg_z[j][k] * vector_mg_z[j][k]);
            }
        }
    }

    // 计算余弦的平均值(MG)
    for (i = 1; i <= NP - 1; i++)
    {
        for (j = 1; j <= NP - 1; j++)
        {
            for (k = 1; k <= t; k++)
            {
                avg_ang_mg[i][j] += ang_mg[i][j][k];
            }
            avg_ang_mg[i][j] = avg_ang_mg[i][j] / t;
        }
    }

    // 计算余弦的标准差(MG)
    for (i = 1; i <= NP - 1; i++)
    {
        for (j = 1; j <= NP - 1; j++)
        {
            for (k = 1; k <= t; k++)
            {
                std_ang_mg[i][j] += (ang_mg[i][j][k] - avg_ang_mg[i][j]) * (ang_mg[i][j][k] - avg_ang_mg[i][j]);
            }
            std_ang_mg[i][j] = sqrt(std_ang_mg[i][j] / t);
        }
    }

    // 打印结果
    for (i = 1; i <= NP - 1; i++)
        {
            for (j = 1; j <= NP - 1; j++)
            {        
                outfile3 << i << " " << j << " " << std_ang_k[i][j]<< endl;
                outfile4 << i << " " << j << " " << std_ang_mg[i][j]<< endl;
            }
        outfile3 << endl;
        outfile4 << endl;
        }
    outfile3.close();
    outfile4.close();

    return 0;
}
