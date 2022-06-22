#include<iostream>
#include<fstream>
#include<cmath>
#define N0 70
#define flag 5000
using namespace std;

float disp[60000][80];

struct example_info
{
    char name[5];
    float x, y, z;
};

example_info p_info[60000][80];

int main()
{
    ifstream infile;
    infile.open("p_60000.pdb");

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
        infile >> p_info[j][k].name >> p_info[j][k].x >> p_info[j][k].y >> p_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    ofstream outfile, outfile1, outfile2;
    outfile.open("trans_time_series.dat", ios::out | ios::ate);
    outfile1.open("trans_time_residue.dat", ios::out | ios::ate);
    outfile2.open("trans_time_example.dat", ios::out | ios::ate);

    // 计算每个磷原子位移随时间步长的变化
    for (i = 1; i <= flag; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            for (k = 1; k <= t - i; k++)
            {
                disp[i][j] += sqrt((p_info[k][j].x - p_info[k + i][j].x) * (p_info[k][j].x - p_info[k + i][j].x) + (p_info[k][j].y - p_info[k + i][j].y) * (p_info[k][j].y - p_info[k + i][j].y) + (p_info[k][j].z - p_info[k + i][j].z) * (p_info[k][j].z - p_info[k + i][j].z));
            }
            disp[i][j] /= (t - i);
            outfile << i << " " << j << " " << disp[i][j] << endl;
        }
        outfile << endl;
    }
    outfile.close();

    float avg_series, temp;
    for (i = 1; i <= N0; i++)
    {
        for (j = 1; j <= flag; j++)
        {
            avg_series = 0;
            for (k = j + 1; k <= flag; k++)
            {
                avg_series += disp[k][i];
            }
            avg_series /= flag - j;
            temp = (disp[j][i] - avg_series) / avg_series;
            if (temp > - 0.01 && temp < 0.01)
            {
                outfile1 << i << " " << j * 0.02 << endl;
                break;
            }            
        }
    }
    outfile1.close();

    for (i = 1; i <= flag; i++)
    {
        outfile2 << i << " " << disp[i][15] << endl;
    }
    outfile2.close();

    return 0;
}
