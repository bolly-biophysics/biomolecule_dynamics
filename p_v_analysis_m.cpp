#include<iostream>
#include<fstream>
#include<cmath>
#define N0 70
using namespace std;

float v_mag[70000][80], v_mag_avg[80], v_mag_std[80];
float v_ang[70000][80], v_ang_avg[80], v_ang_std[80];
float v_ang_ref[70000][80], v_ang_ref_avg[80], v_ang_ref_std[80];
float cos_v[80][80][70000], cos_v_avg[80][80], cos_v_std[80][80];
float mag_corr[80][80], ang_corr[80][80], dynamic_corr[80][80];
float x_avg[80], y_avg[80], z_avg[80];

struct example_info
{
    char name[5];
    float x, y, z;
};

example_info p_info[70000][80];

struct example_v
{
    float x, y, z;
};

example_v v[70000][80];

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

    ofstream outfile_1, outfile_2, outfile_3, outfile_4, outfile_5, outfile_6, outfile_7, outfile_8, outfile_9, outfile_10, outfile_11;
    outfile_1.open("pp_v_cos.dat", ios::out | ios::ate);
    outfile_2.open("p_v_mag_series.dat", ios::out | ios::ate);
    outfile_3.open("p_v_mag_avg_std.dat", ios::out | ios::ate);
    outfile_4.open("pp_v_mag_corr.dat", ios::out | ios::ate);
    outfile_5.open("p_v_ang_series.dat", ios::out | ios::ate);
    outfile_6.open("p_v_ang_avg_std.dat", ios::out | ios::ate);
    outfile_7.open("pp_v_ang_corr.dat", ios::out | ios::ate);
    outfile_8.open("p_v_ang_ref_series.dat", ios::out | ios::ate);
    outfile_9.open("p_v_ang_ref_avg_std.dat", ios::out | ios::ate);
    outfile_10.open("p_v_ang_ref_fb_freq.dat", ios::out | ios::ate);
    outfile_11.open("pp_dynamic_corr.dat", ios::out | ios::ate);

    // 求分子平均结构
    for (i = 1; i <= N0; i++)
    {
        for (j = 1; j <= t; j++)
        {
            x_avg[i] += p_info[j][i].x;
            y_avg[i] += p_info[j][i].y;
            z_avg[i] += p_info[j][i].z;
        }
        x_avg[i] = x_avg[i] / t;
        y_avg[i] = y_avg[i] / t;
        z_avg[i] = z_avg[i] / t;
    }

    // 计算单个磷原子运动速度矢量
    for (i = 1; i <= t - 1; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            v[i][j].x = p_info[i + 1][j].x - p_info[i][j].x;
            v[i][j].y = p_info[i + 1][j].y - p_info[i][j].y;
            v[i][j].z = p_info[i + 1][j].z - p_info[i][j].z;
        }
    }

    // 计算两两磷原子之间运动速度夹角方向随时间的变化
    for (i = 1; i <= N0; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            for (k = 1; k <= t - 1; k++)
            {
                cos_v[i][j][k] = (v[k][i].x * v[k][j].x + v[k][i].y * v[k][j].y + v[k][i].z * v[k][j].z) / (sqrt(v[k][i].x * v[k][i].x + v[k][i].y * v[k][i].y + v[k][i].z * v[k][i].z) * sqrt(v[k][j].x * v[k][j].x + v[k][j].y * v[k][j].y + v[k][j].z * v[k][j].z));
            }
        }
    }

    //计算两两磷原子之间运动速度夹角方向的平均值(代表方向关联性)
    for (i = 1; i <= N0; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            cos_v_avg[i][j] = 0;
            for (k = 1; k <= t - 1; k++)
            {
                cos_v_avg[i][j] += cos_v[i][j][k];
            }
            cos_v_avg[i][j] /= t - 1;
        }
    }

    //计算两两磷原子之间运动速度夹角方向的标准差
    for (i = 1; i <= N0; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            cos_v_std[i][j] = 0;
            for (k = 1; k <= t - 1; k++)
            {
                cos_v_std[i][j] += (cos_v[i][j][k] - cos_v_avg[i][j]) * (cos_v[i][j][k] - cos_v_avg[i][j]);
            }
            cos_v_std[i][j] = sqrt(cos_v_std[i][j] / (t - 1));
            outfile_1 << i << " " << j << " " << cos_v_avg[i][j] << " " << cos_v_std[i][j] << endl;
        }
        outfile_1 << endl;
    }
    outfile_1.close();

    // 计算单个磷原子运动速度大小并输出其随时间的变化
    for (i = 1; i <= t - 1; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            v_mag[i][j] = sqrt(v[i][j].x * v[i][j].x + v[i][j].y * v[i][j].y + v[i][j].z * v[i][j].z);
            outfile_2 << i << " " << j << " " << v_mag[i][j] << endl;
        }
        outfile_2 << endl;
    }

    // 计算单个磷原子运动速度大小的平均值
    for (i = 1; i <= N0; i++)
    {
        v_mag_avg[i] = 0;
        for (j = 1; j <= t - 1; j++)
        {
            v_mag_avg[i] += v_mag[j][i];
        }
        v_mag_avg[i] /= t - 1;
    }

    // 计算单个磷原子运动速度大小的标准差
    for (i = 1; i <= N0; i++)
    {
        v_mag_std[i] = 0;
        for (j = 1; j <= t - 1; j++)
        {
            v_mag_std[i] += (v_mag[j][i] - v_mag_avg[i]) * (v_mag[j][i] - v_mag_avg[i]);
        }
        v_mag_std[i] = sqrt(v_mag_std[i] / (t - 1));
    }

    for (i = 1; i <= N0; i++)
    {
        outfile_3 << i << " " << v_mag_avg[i] << " " << v_mag_std[i] << endl;
    }
    outfile_3.close();

    // 计算两两磷原子之间运动速度的大小关联性
    for (i = 1; i <= N0; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            mag_corr[i][j] = 0;
            for (k = 1; k <= t - 1; k++)
            {
                mag_corr[i][j] += (v_mag[k][i] - v_mag_avg[i]) * (v_mag[k][j] - v_mag_avg[j]);
            }
            mag_corr[i][j] = mag_corr[i][j] / (t - 1) / v_mag_std[i] / v_mag_std[j];
            outfile_4 << i << " " << j << " " << mag_corr[i][j] << endl;
        }
        outfile_4 << endl;
    }
    outfile_4.close();

    // 计算单个磷原子运动速度方向随时间的变化
    for (i = 1; i <= t - 2; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            v_ang[i][j] = (v[i][j].x * v[i + 1][j].x + v[i][j].y * v[i + 1][j].y + v[i][j].z * v[i + 1][j].z) / v_mag[i][j] / v_mag[i + 1][j];
            outfile_5 << i << " " << j << " " << v_ang[i][j] << endl;
        }
        outfile_5 << endl;
    }
    outfile_5.close();

    // 计算单个磷原子运动速度方向的平均值
    for (i = 1; i <= N0; i++)
    {
        v_ang_avg[i] = 0;
        for (j = 1; j <= t - 2; j++)
        {
            v_ang_avg[i] += v_ang[j][i];
        }
        v_ang_avg[i] /= t - 2;
    }

    // 计算单个磷原子运动速度方向的标准差
    for (i = 1; i <= N0; i++)
    {
        v_ang_std[i] = 0;
        for (j = 1; j <= t - 2; j++)
        {
            v_ang_std[i] += (v_ang[j][i] - v_ang_avg[i]) * (v_ang[j][i] - v_ang_avg[i]);
        }
        v_ang_std[i] = sqrt(v_ang_std[i] / (t - 2));
    }

    for (i = 1; i <= N0; i++)
    {
        outfile_6 << i << " " << v_ang_avg[i] << " " << v_ang_std[i] << endl;
    }
    outfile_6.close();

    // 计算两两磷原子之间运动速度方向大小的关联性
    for (i = 1; i <= N0; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            ang_corr[i][j] = 0;
            for (k = 1; k <= t - 2; k++)
            {
                ang_corr[i][j] += (v_ang[k][i] - v_ang_avg[i]) * (v_ang[k][j] - v_ang_avg[j]);
            }
            ang_corr[i][j] = ang_corr[i][j] / (t - 2) / v_ang_std[i] / v_ang_std[j];
            outfile_7 << i << " " << j << " " << ang_corr[i][j] << endl;
        }
        outfile_7 << endl;
    }
    outfile_7.close();

    // 计算单个磷原子运动速度方向随时间的变化(相对平均位置)
    for (i = 1; i <= t - 1; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            float ref_x = p_info[i][j].x - x_avg[j];
            float ref_y = p_info[i][j].y - y_avg[j];
            float ref_z = p_info[i][j].z - z_avg[j];
            v_ang_ref[i][j] = (v[i][j].x * ref_x + v[i][j].y * ref_y + v[i][j].z * ref_z) / v_mag[i][j] / sqrt(ref_x * ref_x + ref_y * ref_y + ref_z * ref_z);
            outfile_8 << i << " " << j << " " << v_ang_ref[i][j] << endl;
        }
        outfile_8 << endl;
    }
    outfile_8.close();

    // 计算单个磷原子运动速度方向的平均值(相对平均位置)
    for (i = 1; i <= N0; i++)
    {
        v_ang_ref_avg[i] = 0;
        for (j = 1; j <= t - 1; j++)
        {
            v_ang_ref_avg[i] += v_ang_ref[j][i];
        }
        v_ang_ref_avg[i] /= t - 1;
    }

    // 计算单个磷原子运动速度方向的标准差(相对平均位置)
    for (i = 1; i <= N0; i++)
    {
        v_ang_ref_std[i] = 0;
        for (j = 1; j <= t - 1; j++)
        {
            v_ang_ref_std[i] += (v_ang_ref[j][i] - v_ang_ref_avg[i]) * (v_ang_ref[j][i] - v_ang_ref_avg[i]);
        }
        v_ang_ref_std[i] = sqrt(v_ang_ref_std[i] / (t - 1));
    }

    for (i = 1; i <= N0; i++)
    {
        outfile_9 << i << " " << v_ang_ref_avg[i] << " " << v_ang_ref_std[i] << endl;
    }
    outfile_9.close();

    // 统计磷原子运动速度(相对平均位置)正向或反向的频率
    float forward[80], backward[80];
    for (i = 1; i <= N0; i++)
    {
        forward[i] = 0, backward[i] = 0;
        for (j = 1; j <= t - 1; j++)
        {
            if (v_ang_ref[j][i] > 0.8)
                forward[i]++;
            else if (v_ang_ref[j][i] < -0.8)
                backward[i]++;
        }
        outfile_10 << i << " " << forward[i] / (t - 1) << " " << backward[i] / (t - 1) << endl;
    }
    outfile_10.close();

    // 计算两两磷原子之间动态关联特性
    float avg_ref_i, avg_ref_j;
    for (i = 1; i <= N0; i++)
    {
        for (j = 1; j <= N0; j++)
        {
            dynamic_corr[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                float ref_i_x = p_info[k][i].x - x_avg[i];
                float ref_i_y = p_info[k][i].y - y_avg[i];
                float ref_i_z = p_info[k][i].z - z_avg[i];
                float ref_j_x = p_info[k][j].x - x_avg[j];
                float ref_j_y = p_info[k][j].y - y_avg[j];
                float ref_j_z = p_info[k][j].z - z_avg[j];
                dynamic_corr[i][j] += (ref_i_x * ref_j_x + ref_i_y * ref_j_y + ref_i_z * ref_j_z);
                avg_ref_i += (ref_i_x * ref_i_x + ref_i_y * ref_i_y + ref_i_z * ref_i_z);
                avg_ref_j += (ref_j_x * ref_j_x + ref_j_y * ref_j_y + ref_j_z * ref_j_z);
            }
            avg_ref_i = sqrt(avg_ref_i / t);
            avg_ref_j = sqrt(avg_ref_j / t);
            dynamic_corr[i][j] = dynamic_corr[i][j] / t / avg_ref_i / avg_ref_j;
            outfile_11 << i << " " << j << " " << dynamic_corr[i][j] << endl;
        }
        outfile_11 << endl;
    }
    outfile_11.close();

    return 0;
}
