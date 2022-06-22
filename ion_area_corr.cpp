#include<iostream>
#include<fstream>
#include<cmath>
#define NP 70
#define NMG 104
#define R 6
using namespace std;

float n_mg[60000][80], avg_mg[80], std_mg[80], corr[80][80];

struct example_info
{
    char name[5];
    float x, y, z;
};

example_info p_info[60000][80];
example_info mg_info[60000][120];

int main()
{
    ifstream infile1, infile2;
    infile1.open("p.pdb");
    infile2.open("mg.pdb");

    if (!infile1)
    {
        cerr << "Open infile1 failure!" << endl;
        return -1;
    }

    // 将磷原子坐标格式化存入结构数组
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
        infile1 >> p_info[j][k].name >> p_info[j][k].x >> p_info[j][k].y >> p_info[j][k].z;
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

    // 将镁离子坐标格式化存入结构数组
    i = 1;
    while (!infile2.eof())
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
        infile2 >> mg_info[j][k].name >> mg_info[j][k].x >> mg_info[j][k].y >> mg_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile2.close();

    ofstream outfile1, outfile2;
    outfile1.open("mg_avg_std.dat", ios::out | ios::ate);
    outfile2.open("ion_area_corr.dat", ios::out | ios::ate);

    // 计算每一时刻单个磷原子附近的镁离子个数
    int frame;
    float d;
    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= NP; j++)
        {
            frame = 0;
            for (k = 1; k <= NMG; k++)
            {
                d = sqrt((mg_info[i][k].x - p_info[i][j].x) * (mg_info[i][k].x - p_info[i][j].x) + (mg_info[i][k].y - p_info[i][j].y) * (mg_info[i][k].y - p_info[i][j].y) + (mg_info[i][k].z - p_info[i][j].z) * (mg_info[i][k].z - p_info[i][j].z));
                if (d < R)
                    frame++;
            }
            n_mg[i][j] = frame;
        }
    }

    // 计算单个磷原子附近的平均镁离子个数
    for (i = 1; i <= NP; i++)
    {
        avg_mg[i] = 0;
        for (j = 1; j <= t; j++)
        {
            avg_mg[i] += n_mg[j][i];
        }
        avg_mg[i] = avg_mg[i] / t;
    }

    // 计算单个磷原子附近镁离子个数涨落的标准差
    for (i = 1; i <= NP; i++)
    {
        std_mg[i] = 0;
        for (j = 1; j <= t; j++)
        {
            std_mg[i] += (n_mg[j][i] - avg_mg[i]) * (n_mg[j][i] - avg_mg[i]);
        }
        std_mg[i] = sqrt(std_mg[i] / t);
    }

    for (i = 1; i <= NP; i++)
    {
        outfile1 << i << " " << avg_mg[i] << " " << std_mg[i] << endl;
    }
    outfile1.close();

    // 计算两两磷原子附近镁离子个数之间的相关系数
    for (i = 1; i <= NP; i++)
        {
            for (j = 1; j <= NP; j++)
            {        
                corr[i][j] = 0;
                for (k = 1; k <= t; k++)
                {
                    corr[i][j] += (n_mg[k][i] - avg_mg[i]) * (n_mg[k][j] - avg_mg[j]);
                }
            corr[i][j] = corr[i][j] / t / std_mg[i] / std_mg[j];
            }
        }   

    // 打印结果
    for (i = 1; i <= NP; i++)
        {
            for (j = 1; j <= NP; j++)
            {        
                outfile2 << i << " " << j << " " << corr[i][j] << endl;
            }
        outfile2 << endl;
        }
    outfile2.close();

    return 0;
}
