#include<iostream>
#include<fstream>
#include<cmath>
#define p1_len 17
#define p3_len 11
using namespace std;

float axis_angle[7000], d[7000];

struct axis_info
{
    float x, y, z;
};

axis_info p1_info[7000][20];
axis_info p3_info[7000][20];

struct axis_vector
{
    float x, y, z;
};

axis_vector p1_vector[7000];
axis_vector p3_vector[7000];

int main()
{
    ifstream infile;
    infile.open("axis_p1.dat");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将中心轴坐标格式化存入结构数组
    int i = 1, j, k, t;
    float aa;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (p1_len + 0.0); 
        k = fmod(i, p1_len);
        if (k == 0) 
        {
            k = p1_len;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> p1_info[j][k].x >> p1_info[j][k].y >> p1_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    infile.open("axis_p3.dat");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    i = 1;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (p3_len + 0.0); 
        k = fmod(i, p3_len);
        if (k == 0) 
        {
            k = p3_len;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> p3_info[j][k].x >> p3_info[j][k].y >> p3_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    ofstream outfile;
    outfile.open("axis_angle.dat", ios::out | ios::ate);

    // 提取p1方向矢量
    for (i = 1; i <= t; i++)
    {
        p1_vector[i].x = p1_info[i][p1_len - 2].x - p1_info[i][3].x;
        p1_vector[i].y = p1_info[i][p1_len - 2].y - p1_info[i][3].y;
        p1_vector[i].z = p1_info[i][p1_len - 2].z - p1_info[i][3].z;
    }

    // 提取p3方向矢量
    for (i = 1; i <= t; i++)
    {
        p3_vector[i].x = p3_info[i][p3_len - 2].x - p3_info[i][3].x;
        p3_vector[i].y = p3_info[i][p3_len - 2].y - p3_info[i][3].y;
        p3_vector[i].z = p3_info[i][p3_len - 2].z - p3_info[i][3].z;
    }

    // 计算p1与p3方向矢量之间的夹角
    axis_angle[i] = 0; d[i] = 0;
    for (i = 1; i <= t; i++)
    {
        axis_angle[i] = (p1_vector[i].x * p3_vector[i].x + p1_vector[i].y * p3_vector[i].y + p1_vector[i].z * p3_vector[i].z) / (sqrt(p1_vector[i].x * p1_vector[i].x + p1_vector[i].y * p1_vector[i].y + p1_vector[i].z * p1_vector[i].z) * sqrt(p3_vector[i].x * p3_vector[i].x + p3_vector[i].y * p3_vector[i].y + p3_vector[i].z * p3_vector[i].z));
        axis_angle[i] = acos(axis_angle[i]) / 3.14 * 180;
        d[i] = sqrt((p1_info[i][p1_len].x - p3_info[i][1].x) * (p1_info[i][p1_len].x - p3_info[i][1].x) + (p1_info[i][p1_len].y - p3_info[i][1].y) * (p1_info[i][p1_len].y - p3_info[i][1].y) + (p1_info[i][p1_len].z - p3_info[i][1].z) * (p1_info[i][p1_len].z - p3_info[i][1].z));
        outfile << i << " " << axis_angle[i] << " " << d[i] << endl;
    }
    outfile.close();

    return 0;
}
