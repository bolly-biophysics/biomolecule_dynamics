#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#define N_RNA_atom 11587
#define N_MG_add 515
using namespace std;

struct Molecule_info
{
    string ATOM, atomName, resName, chainID, element;
    int serial, resSeq;
    float x, y, z, ocpy, bfac;
};

Molecule_info RNA[10000];

int main()
{
    ifstream infile1;
    infile1.open("1y26_pure.pdb", ios::in);

    if (!infile1)
    {
        cerr << "Open infile1 failure!" << endl;
        return -1;
    }

    // 将适体信息格式化存入结构数组
    int i = 1, j, k, t;
    while (!infile1.eof())
    {
        infile1 >> RNA[i].ATOM >> RNA[i].serial >> RNA[i].atomName >> RNA[i].resName >> RNA[i].chainID >> RNA[i].resSeq >> RNA[i].x >> RNA[i].y >> RNA[i].z >> RNA[i].ocpy >> RNA[i].bfac >> RNA[i].element;
        i++;
    }
    t = i - 2;
    cout << t << endl;
    infile1.close();

    // 计算适体几何中心坐标
    float RNA_center_x, RNA_center_y, RNA_center_z;
    for (i = 1; i <= N_RNA_atom; i++)
    {
        RNA_center_x += RNA[i].x;
        RNA_center_y += RNA[i].y;
        RNA_center_z += RNA[i].z;
    }
    RNA_center_x /= N_RNA_atom;
    RNA_center_y /= N_RNA_atom;
    RNA_center_z /= N_RNA_atom;

    // 计算适体边界最大距离
    float d_max = 10;
    for (i = 1; i <= N_RNA_atom; i++)
    {
        float d = sqrt((RNA[i].x - RNA_center_x) * (RNA[i].x - RNA_center_x) + (RNA[i].y - RNA_center_y) * (RNA[i].y - RNA_center_y) + (RNA[i].z - RNA_center_z) * (RNA[i].z - RNA_center_z));
        if (d > d_max)
            d_max = d;
    }

    // 原适体PDB文件末尾添加TER
    ofstream outfile;
    outfile.open("1y26_MG.pdb", ios::out | ios::app);
    outfile << "TER" << endl;
    outfile.close();

    // 在给定的限制条件下随机产生新配体的几何中心坐标，并将其完整信息追加写入适体的PDB文件
    float Lig_center_tmp_x, Lig_center_tmp_y, Lig_center_tmp_z;
    float d_Lig_RNA, d_ij;
    float Lig_center_loc_x[1000], Lig_center_loc_y[1000], Lig_center_loc_z[1000];
    srand((unsigned int)time(NULL));
    for (i = 1; i <= N_Lig_add; i++)
    {
        bool gen_center_failed = true;
        while (gen_center_failed)
        {
            bool is_too_close = false;
            Lig_center_tmp_x = rand() % 2000 - 1000;
            Lig_center_tmp_y = rand() % 2000 - 1000;
            Lig_center_tmp_z = rand() % 2000 - 1000;
            d_Lig_RNA = sqrt((Lig_center_tmp_x - RNA_center_x) * (Lig_center_tmp_x - RNA_center_x) + (Lig_center_tmp_y - RNA_center_y) * (Lig_center_tmp_y - RNA_center_y) + (Lig_center_tmp_z - RNA_center_z) * (Lig_center_tmp_z - RNA_center_z));
            if (d_Lig_RNA > d_max + 5 && d_Lig_RNA < d_max + 10)
            {
                Lig_center_loc_x[i] = Lig_center_tmp_x;
                Lig_center_loc_y[i] = Lig_center_tmp_y;
                Lig_center_loc_z[i] = Lig_center_tmp_z;
            }
            else
            {
                continue;
            }
            for (j = 1; j < i; j++)
            {
                d_ij = sqrt((Lig_center_loc_x[i] - Lig_center_loc_x[j]) * (Lig_center_loc_x[i] - Lig_center_loc_x[j]) + (Lig_center_loc_y[i] - Lig_center_loc_y[j]) * (Lig_center_loc_y[i] - Lig_center_loc_y[j]) + (Lig_center_loc_z[i] - Lig_center_loc_z[j]) * (Lig_center_loc_z[i] - Lig_center_loc_z[j]));
                if (d_ij < 5)
                {
                    is_too_close = true;
                    break;
                }
            }
            if (is_too_close)
            {
                continue;
            }
            else
            {
                outfile.open("1y26_MG.pdb", ios::out | ios::app);
                outfile << setw(6) << "HETATM" << setw(5) << i << "  " << setw(3) << setiosflags(ios::left) << "MG" << " " << setw(3) << resetiosflags(ios::left) << "MG" << " " << "X" << setw(4) << resetiosflags(ios::right) << i << "    " << setw(8) << fixed << setprecision(3) << Lig_center_loc_x[i] << setw(8) << fixed << setprecision(3) << Lig_center_loc_y[i] << setw(8) << fixed << setprecision(3) << Lig_center_loc_z[i] << setw(6) << fixed << setprecision(2) << "1.00" << setw(6) << fixed << setprecision(2) << "0.00" << "          " << setw(2) << "MG" << endl;
                outfile.close();
            }
            gen_center_failed = false;
        }
    }

    // 适配体PDB文件末尾添加END
    outfile.open("1y26_MG.pdb", ios::out | ios::app);
    outfile << "END" << endl;
    outfile.close();

    return 0;
}
