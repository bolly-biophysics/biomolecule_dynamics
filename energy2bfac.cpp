#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
using namespace std;

struct example_info
{
    string ATOM, atomType, resName, element;
    int serial, resSeq;
    float x, y, z, ocpy, bfac;
};

example_info pdb[3000];

int main()
{
    ifstream infile, infile1;
    infile.open("50.pdb");
    infile1.open("energy.dat");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将PDB信息格式化存入结构数组
    int i = 1, t;
    while (!infile.eof())
    {
        infile >> pdb[i].ATOM >> pdb[i].serial >> pdb[i].atomType >> pdb[i].resName >> pdb[i].resSeq >> pdb[i].x >> pdb[i].y >> pdb[i].z >> pdb[i].ocpy >> pdb[i].bfac >> pdb[i].element;
        i++;
    }
    t = i - 2;
    cout << t << endl;
    infile.close();

    if (!infile1)
    {
        cerr << "Open infile1 failure!" << endl;
        return -1;
    }

    // 将各残基对结合自由能的贡献存入数组
    i = 1;
    float energy[80];
    while (!infile1.eof())
    {
        infile1 >> energy[i];
        i++;
    }
    int n = i - 2;
    cout << n << endl;
    infile1.close();

    ofstream outfile;
    outfile.open("energy.pdb", ios::out | ios::ate);

    for (i = 1; i <= t; i++)
    {
        if (pdb[i].atomType.length() <= 3)
        outfile << setw(4) << pdb[i].ATOM << "  " << setw(5) << pdb[i].serial << "  " << setiosflags(ios::left) << setw(3) << pdb[i].atomType << resetiosflags(ios::left) << " " << setw(3) << pdb[i].resName << "  " << setw(4) << pdb[i].resSeq << "    " << setw(8) << pdb[i].x << setw(8) << pdb[i].y << setw(8) << pdb[i].z << setw(6) << "1.00" << setw(6) << energy[pdb[i].resSeq] << "          " << setw(2) << pdb[i].element << endl;
        else
        outfile << setw(4) << pdb[i].ATOM << "  " << setw(5) << pdb[i].serial << " " << setw(4) << pdb[i].atomType << " " << setw(3) << pdb[i].resName << "  " << setw(4) << pdb[i].resSeq << "    " << setw(8) << pdb[i].x << setw(8) << pdb[i].y << setw(8) << pdb[i].z << setw(6) << "1.00" << setw(6) << energy[pdb[i].resSeq] << "          " << setw(2) << pdb[i].element << endl;
    }
    outfile.close();

    return 0;
}
