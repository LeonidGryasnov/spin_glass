#include <iostream>
#include <fstream>
#include <string>
constexpr int count = 3000;
constexpr int averaging = 100;
double data[count][2];

using namespace std;

string get_file_name(int i, string s)
{
    return to_string(i) + s;
}

void init_data(double data_arr[count][2])
{
    for(auto i = 0; i < count; i++){
        for(auto j = 0; j < 2; j++){
             data_arr[i][j] = 0;
        }
}
}

int read_file(int file_number)
{
     string s = " - Xsg(t,L).dat";
    double arr[count][2];
    string file_name = get_file_name(file_number, s);
    ifstream file(file_name);

    if (!file.is_open())
    {
        cout << "Файл не может быть окрыт или создан!" << endl;
        return -1;
    }

    init_data(arr);

    for(auto i = 0; i < count; i++){
        for(auto j = 0; j < 2; j++){
            file >> arr[i][j];
            data[i][j] += arr[i][j];
        }
    }
    return 0;

}

void averaging_data(double arr[count][2])
{
    for(auto i = 0; i < count; i++){
        for(auto j = 0; j < 2; j++){
            arr[i][j] /= averaging;
        }
    }
}

int write_file(string s)
{
    ofstream file(s);

    if (!file.is_open())
    {
        cout << "Файл не может быть окрыт или создан!" << endl;
        return -1;
    }

     for(auto i = 0; i < count; i++){
        file << '\n';
        for(auto j = 0; j < 2; j++){
            file << data[i][j] << '\t';
        }
     }
    return 0;
}

int main()
{
    string s = "L=6;_Stat_Xsg.dat";

    init_data(data);

    for(auto i = 0; i < averaging; i++){
        read_file(i);
    }

    averaging_data(data);
    write_file(s);

    for(auto i = 0; i < count; i++){
        cout << '\n';
        for(auto j = 0; j < 2; j++){
            cout << data[i][j] << '\t';
        }
    }

    return 0;
}
