// gorenie_gauss.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

double Q() {
    return pow(10, 6);
}

double Cp(double T)
{
    return 5 + 0.1 * T + pow(10, -3) * pow(T, 2);
}

double Lambda(double T) {
    return 1 + 0.01 * T + pow(10, -5) * pow(T, 2);
    //return 10;
}

double Lambda_prime(double T) {
    return 0.01 / 2. + pow(10, -5) * T / 2.;
    //return 0;
}

double FprimeLeft(vector<double>& res_vect, vector<double>& x,const int i)
{
    double h_left = x[i] - x[i - 1];
    double h = x[i + 1] - x[i];
    return -(2. / (x[i + 1] - x[i - 1]))* (-Lambda_prime((res_vect[i] + res_vect[i - 1]) / 2.) * (res_vect[i] - res_vect[i - 1]) / (x[i] - x[i - 1])
        + Lambda((res_vect[i] + res_vect[i - 1]) / 2.) / (x[i] - x[i - 1])) 
        - h / h_left / (h + h_left) * res_vect[res_vect.size() - 1] * Cp(res_vect[i]);
}

double FprimeCenter(vector<double>& res_vect, vector<double>& x, const int i)
{
    double h_left = x[i] - x[i - 1];
    double h = x[i + 1] - x[i];
    return -(2. / (x[i + 1] - x[i - 1])) * (Lambda_prime((res_vect[i] + res_vect[i + 1]) / 2.) * (res_vect[i + 1] - res_vect[i]) / (x[i + 1] - x[i])
        - Lambda((res_vect[i] + res_vect[i + 1]) / 2.) / (x[i + 1] - x[i]) - Lambda_prime((res_vect[i] + res_vect[i - 1]) / 2.) * (res_vect[i] - res_vect[i - 1]) / (x[i] - x[i - 1])
        - Lambda((res_vect[i] + res_vect[i - 1]) / 2.) / (x[i] - x[i - 1]))
       + (h - h_left) / h_left / h * res_vect[res_vect.size() - 1] * Cp(res_vect[i]);
}

double FprimeRight(vector<double>& res_vect, vector<double>& x, const int i)
{
    double h_left = x[i] - x[i - 1];
    double h = x[i + 1] - x[i];
    //cout << "m = " << res_vect[res_vect.size() - 1] << endl;
    return -(2. / (x[i + 1] - x[i - 1])) * (Lambda_prime((res_vect[i + 1] + res_vect[i]) / 2.) * (res_vect[i + 1] - res_vect[i]) / (x[i + 1] - x[i])
        + Lambda((res_vect[i + 1] + res_vect[i]) / 2.) / (x[i + 1] - x[i]))
        + h_left / h / (h + h_left) * res_vect[res_vect.size() - 1] * Cp(res_vect[i]);
   
}

double FprimeM(vector<double>& res_vect, vector<double>& x, const int i)
{
    double h_left = x[i] - x[i - 1];
    double h = x[i + 1] - x[i];
    return (h_left / (h * (h + h_left)) * res_vect[i + 1] + (h - h_left) / h / h_left * res_vect[i]
        - h / h_left / (h + h_left) * res_vect[i - 1]) * Cp(res_vect[i]);
}

void F_vector(vector<double>& f, vector<double>& T, const int N1, const int N_M, vector<double>& x,int N_start, int N_finish) {
    double h_left;
    double h;
    double q_tmp = 0;
    for (int j = 1; j <= N1; j++)
    {
        if (j >= N_start && j <= N_finish)
        {
            q_tmp = Q();
            //cout << "q_tmp = " << q_tmp << endl;

        }
        else {
            q_tmp = 0.;
        }
        h_left = x[j] - x[j - 1];
        h = x[j + 1] - x[j];
        f[j - 1] = (2. / (x[j + 1] - x[j - 1])) * (Lambda((T[j + 1] + T[j]) / 2.) * (T[j + 1] - T[j]) / (x[j + 1] - x[j]) - Lambda((T[j] + T[j - 1]) / 2.) * (T[j] - T[j - 1]) / (x[j] - x[j - 1]))
            - Cp(T[j]) * T[N_M] * (h_left / h / (h + h_left) * T[j + 1] + (h - h_left) / h / h_left * T[j] - h / h_left / (h + h_left) * T[j - 1]) 
            + q_tmp;
    }
}

void Jacobian(vector<vector<double>>& jac, vector<double>& T, const int N1,const int N_M, vector<double>& x)
{
    // уменьшили количество уравнений
    T[1] = T[0];
    T[N_M - 1] = T[N_M - 2]; //адиабатическое граничное условие
    jac[0][0] = FprimeRight(T, x, 1); // dF2/dT3
    jac[0][N1 - 1] = FprimeM(T, x, 1); //dF2/dM
    jac[1][0] = FprimeCenter(T, x, 2); //dF3/dT3
    jac[1][1] = FprimeRight(T, x, 2); //dF3/dT4
    jac[1][N1 - 1] = FprimeM(T, x, 2); //dF3/dM
    for (int i = 2; i < N1 - 1; i++) {
        jac[i][i - 2] = FprimeLeft(T, x, i + 1); 
        jac[i][i - 1] = FprimeCenter(T, x, i + 1);
        jac[i][i] = FprimeRight(T, x, i + 1);
        jac[i][N1 - 1] = FprimeM(T, x, i + 1);
    }
    jac[N1 - 1][N1 - 1 - 2] = FprimeLeft(T, x, N1 - 1 + 1);
    jac[N1 - 1][N1 - 1 - 1] = FprimeCenter(T, x, N1 - 1 + 1) + FprimeRight(T, x, N1 - 1 + 1);
    jac[N1 - 1][N1 - 1] = FprimeM(T, x, N1 - 1 + 1);
}

void gauss(vector<vector<double>>& a, vector<double>& y, int n, vector<double>& x)
{
    double max;
    int k, index;
    const double eps = 0.000001; 
    k = 0;
    while (k < n)
    {
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(a[i][k]) > max)
            {
                max = abs(a[i][k]);
                index = i;
            }
        }
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (abs(temp) < eps) continue; 
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)  continue;
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
        k++;
    }
    for (k = n - 1; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
}

void Solve(const int N, ofstream& fout, vector<double>& x, vector <double>& res_vect)
{
    ofstream file_h;
    file_h.open("file_h" + std::to_string(N) + ".dat");
    const int N_syst = N - 2;
    const int N_M = N;
    double N_start = 0, N_finish = N - 1;
    double eps = (x[N_M - 2] - x[0]) / (N_M - 2.) / 2.;
    double y;
    vector<double> dRes(N_syst), f_nevaz(N_syst);
    vector < vector <double> > jac(N_syst);
    for (int i = 0; i < N_syst; i++) {
        jac[i].resize(N_syst);
    }
    double f_mod = 100;
    double f_mod_tmp = 0;
    int num_iter = 0;
    file_h << "TITLE=\"" << "Graphics" << "\"" << endl;
    file_h << R"(VARIABLES= "i", "F")" << endl;

    for (int j = 0; j < N_M - 1; j++)
    {
        if (fabs(x[j]) > 0.50 && fabs(x[j] - eps) < 0.50) {
            N_start = j;
            cout << "N_start = " << N_start << endl;
        }
        if (fabs(x[j] + eps) > 0.51 && fabs(x[j]) < 0.51) {
            N_finish = j;
            cout << "N_finish = " << N_finish << endl;
        }
    }

    cout << "Number = " << N_finish - N_start + 1 << endl;
    F_vector(f_nevaz, res_vect, N_syst, N_M, x, N_start, N_finish);

    while (fabs(f_mod - f_mod_tmp) > eps)
    {
        f_mod_tmp = f_mod;
        f_mod = 0;
        Jacobian(jac, res_vect, N_syst, N_M, x);
        cout << "iteration N " << num_iter << endl;
        gauss(jac, f_nevaz, N_syst, dRes);
        for (int i = 2; i <= N_M - 2; i++) {
            res_vect[i] += dRes[i - 2];
        }
        res_vect[N_M - 1] = res_vect[N_M - 2];
        cout << "T[" << N_M << "] = " << res_vect[N_M] << endl;
        res_vect[N_M] += dRes[N_syst - 1];
        num_iter += 1;
        F_vector(f_nevaz, res_vect, N_syst, N_M, x, N_start, N_finish);

        for (int i = 0; i < N_syst; i++) {
            f_mod += f_nevaz[i] * f_nevaz[i];
        }
        f_mod = pow(f_mod, 0.5);
        
        
        file_h << num_iter << " " << f_mod << endl;
        cout << endl << "modul F = " << f_mod << " ";
        std::cout << std::endl << endl;;
    }
    
    
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << R"(VARIABLES= "x", "T", "lambda")" << endl;
    
    for (int i = 0; i < N_M; i++) {
        fout << x[i] << " " << res_vect[i] << " " << Lambda(res_vect[i]) << endl;
    }
    std::cout << "Hello World!\n";
    file_h.close();
}

void Init(const int N, double x_l, double x_r, const double T_l, const double T_r, vector<double>& x, vector<double>& T_vect)
{
    double M = 0.5;
    T_vect.resize(N + 1);
    x.resize(N);
    x[0] = x_l;
    //double h = (x_r - x_l) * 2. / (N - 1) / (N);
    double h = (x_r - x_l) / (N - 1);
    T_vect[0] = T_l;
    T_vect[1] = T_l;
    x[1] = h;
    for (int i = 2; i < N; i++)
    {
        T_vect[i] = T_l;
        //x[i] = x[i - 1] + h * i;
        x[i] = h * i;
       // cout << "i = " << i << "  x[i] = " << x[i] << endl;

    }
    T_vect[N - 1] = T_vect[N - 2];
    T_vect[N] = M;
}

int main()
{
    const int N = 2000;
    const double T_l = 500;
    const double T_r = 5000;
    const double x_l = 0;
    const double x_r = 1;
    vector<double> x;
    vector<double> T_vect;
    Init(N, x_l, x_r, T_l, T_r, x, T_vect);
    ofstream fout;
    fout.open("1612file" + std::to_string(N) + ".dat");
    Solve(N, fout, x, T_vect);
    fout.close();
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
