#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <random>
#include <string>
using namespace std;

//генератор рандомной последовательности
vector<char> randomize(vector<char> vec){
    vector<char> res = std::move(vec);
    random_device rd;
    default_random_engine engine(rd());
    shuffle(res.begin(), res.end(), engine);
    return res;
}

//информация между матрицами
double Kullback_information_mat_mat(vector<vector<int>> m1, vector<vector<int>> m2, int alphabet, int n){
    double I = 0;
    for(int j = 0; j < n; j++){
        double I_j = 0;
        for(int i = 0; i < alphabet; i++){
            I_j = I_j + m1[i][j]*log(m1[i][j]);
            I_j = I_j + m2[i][j]*log(m2[i][j]);
            I_j = I_j - (m1[i][j]+m2[i][j])*log(m1[i][j]+m2[i][j]);
        }
        double s1 = 0;
        double s2 = 0;
        for(int k = 0; k < alphabet; k++){
            s1 += m1[k][j];
            s2 += m2[k][j];
        }
        I_j = I_j + ((s1+s2))*log(s1+s2);
        I_j = I_j - s1*log(s1) - s2*log(s2);
        I = I + I_j;
    }
    return I;
}

//мера различия
double measure_of_difference(vector<vector<int>> m1, vector<vector<int>> m2, int alphabet, int n){
    return sqrt(4* Kullback_information_mat_mat(std::move(m1), std::move(m2), alphabet,n)) - sqrt(2*(alphabet-1)*(n-1) - 1);
}

//генерация случайной матрицы из 0 и 1
vector<vector<int>> random_matrix(int alphabet, int n){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);
    vector<vector<int>> matrix(alphabet, vector<int>(n));
    for (int i = 0; i < alphabet; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = dis(gen);
        }
    }
    return matrix;
}

//норма вектора в l2
double norm(const vector<double>& v, int n){
    double sum = 0.0;
    for(double i : v) {
        sum += i * i * n;
    }
    return std::sqrt(sum);
}

//скалярное произведение
double scalar_product(const vector<vector<int>>& m, const vector<double>& prob, int alphabet, int n){
    double q = 0;
    for (int i = 0; i < alphabet; ++i) {
        for (int j = 0; j < n; ++j) {
            q = q+ m[i][j]*prob[i];
        }
    }
    return q;
}

//преобразование случайной матрицы
vector<vector<double>> random_matrix_modified(const vector<vector<int>>& m, int alphabet, int n, vector<double> prob, double A, double B){
    vector<vector<double>> result(alphabet, vector<double>(n));
    vector<vector<double>> m1(alphabet, vector<double>(n));
    vector<vector<double>> m2(alphabet, vector<double>(n));
    for (int i = 0; i < alphabet; i++) {
        for (int j = 0; j < n; j++) {
            double D = (scalar_product(m, prob, alphabet, n) - B)/norm(prob,n);
            m1[i][j] = m[i][j] - (D * prob[i])/norm(prob,n);
            m2[i][j] = (B * prob[i])/pow(norm(prob,n),2);
        }
    }
    double t = 0.0;
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    for (int i = 0; i < alphabet; i++) {
        for (int j = 0; j < n; j++) {
            double k = pow(m1[i][j],2);
            a = a + pow(m2[i][j] - m1[i][j],2);
            b = b + 2 * (m2[i][j]*m1[i][j] - k);
            c = c + k;
        }
    }
    c = c-A;
    double D = pow(b,2) - 4*a*c;
    t = (-b - sqrt(D))/(2*a);
    for (int i = 0; i < alphabet; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = m1[i][j] + t*(m2[i][j] - m1[i][j]);
        }
    }
    return result;
}

//рассчет матрицы F и F^{-1}. В F^{-1} 0 -> 0; 1-> вверх; 2-> влево; 3-> диагональ
pair<vector<vector<double>>, vector<vector<int>>> F(const vector<vector<double>>& m_mod, double d, vector<char> seq, int n){
    int N = (int) seq.size();
    vector<vector<double>> F(N+1, vector<double>(N+1));
    vector<vector<int>> F2(N+1, vector<int>(N+1));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if(i == 0 or j == 0){
                F[i][j] = 0.0;
                F2[i][j] = 0;
            }
            else{
                int k = j - n * int((j-0.1)/n) - 1;
                //cout << "k= "<< k<<" ";
                int s_i = int(seq[i] - 'A');
                double i_j = F[i-1][j] - d;
                double ij_ = F[i][j-1] - d;
                double i_j_ = F[i-1][j-1] +m_mod[s_i][k];
                //cout << "i_j_ = "<< i_j_<<"\n";
                if(abs(i_j) < pow(10,-8) and  abs(ij_) < pow(10,-8) and abs(i_j_) < pow(10,-8)){
                    i_j = 0;
                    ij_ = 0;
                    i_j_ = 0;
                }
                else {
                    if (0.0 > ij_ and 0.0 > i_j and 0.0 > i_j_) {
                        F[i][j] = 0.0;
                        F2[i][j] = 0;
                    } else if (i_j > ij_ and i_j > i_j_) {
                        F[i][j] = i_j;
                        F2[i][j] = 1;
                    } else if (ij_ > i_j_) {
                        F[i][j] = ij_;
                        F2[i][j] = 2;
                    } else {
                        F[i][j] = i_j_;
                        //cout<<"|||"<< s_i<<"  "<< k<< "  "<< m_mod[s_i][k] << "|||"<<"\n";
                        F2[i][j] = 3;
                    }
                }
            }
            //cout<< F2[i][j]<<" ";
        }
    }
    pair<vector<vector<double>>, vector<vector<int>>> result = pair(F,F2);
    return result;
}

//расчет координат максимума F
tuple<int, int, double> maxim(const vector<vector<double>>& F, int N){
    double maximum = 0.0;
    int i_max = 0;
    int j_max = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (F[i][j] > maximum) {
                maximum = F[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }
    tuple<int,int, double> res(i_max, j_max, maximum);
    return res;
}
//разница i_m - i_0
int difference(const pair<vector<vector<double>>, vector<vector<int>>>& pair){
    vector<vector<double>> F = pair.first;
    vector<vector<int>> F2 = pair.second;
    tuple<int, int, double> p1 = maxim(F, (int)F.size());
    int i_m = get<0>(p1);
    int j_m = get<1>(p1);
    int i_0 = i_m;
    int i = i_m;
    int j = j_m;
    while(F2[i][j] != 0){
        if (F2[i][j] != 2){
            i_0 -= 1;
        }
        switch (F2[i][j]) {
            case 1:
                i = i - 1;
                break;
            case 2:
                j = j - 1;
                break;
            case 3:
                j = j - 1;
                i = i - 1;
                break;
        }
    }
    return i_m - i_0;
}

//среднее + дисперсия
double findAverage(const vector<double>& numbers) {
    if (numbers.empty()) {
        throw std::invalid_argument("The list of numbers cannot be empty.");
    }

    double sum = 0.0;
    for (double num : numbers) {
        sum += num;
    }

    return sum / (int)numbers.size();
}
double findVariance(const vector<double>& numbers) {
    if (numbers.empty()) {
        throw std::invalid_argument("The list of numbers cannot be empty.");
    }

    double mean = findAverage(numbers);
    double sumSquaredDifferences = 0.0;

    for (double num : numbers) {
        double difference = num - mean;
        sumSquaredDifferences += difference * difference;
    }

    return sumSquaredDifferences / (int)numbers.size();
}

//генератор последовательности со вставками и делециями
vector<char> generatePeriodicVector(const vector<char>& pattern, int length) {
    vector<char> result;
    int patternLength = (int)pattern.size();
    for (int i = 0; i < length; ++i) {
        result.push_back(pattern[i % patternLength]);
    }

    return result;
}
char getRandomChar(mt19937& rng, int alphabet) {
    std::uniform_int_distribution<int> dist(0, alphabet-1);
    return 'A' + dist(rng);
}
void randomInsertChars(vector<char>& vec, int n, int alphabet) {
    random_device rd;
    mt19937 rng(rd());
    for (int i = 0; i < n; ++i) {
        char randomChar = getRandomChar(rng, alphabet);
        uniform_int_distribution<int> positionDist(0, (int)vec.size());
        int position = positionDist(rng);
        vec.insert(vec.begin() + position, randomChar);
    }
}
void randomRemoveChars(vector<char>& vec, int n) {
    random_device rd;
    mt19937 rng(rd());
    if (n > vec.size()) {
        n = (int)vec.size();
    }
    vector<int> indices;
    uniform_int_distribution<int> dist(0, (int)vec.size() - 1);
    for (int i = 0; i < n; ++i) {
        int randomIndex = dist(rng);
        indices.push_back(randomIndex);
    }
    sort(indices.begin(), indices.end());
    for (int index : indices) {
        vec.erase(vec.begin() + index);
    }
}
vector<char> generator(int len, int number_of_insertion, int number_of_deletion, int alphabet){
    string period = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    vector<char> pattern(period.begin(), period.end());
    vector<char> periodic = generatePeriodicVector(pattern, len);
    randomInsertChars(periodic,number_of_insertion, alphabet);
    randomRemoveChars(periodic,number_of_deletion);
    return periodic;
}

//вектор вероятностей
vector<double> prob(const vector<char>& v, int alphabet){
    vector<double> res(alphabet,0.0);
    for (char i : v) {
        int f = i - 'A';
        res[f] += 1.0/ v.size();
    }
    return res;
}


int main() {
    int alphabet = 26;
    vector<double> res;
    vector<char> s = generator(500,0,0, alphabet);
    vector<char> s11 = randomize(s);
    vector<char> s12 = randomize(s);
    vector<char> sequense;
    sequense.insert(sequense.end(), s11.begin(), s11.end());
    sequense.insert(sequense.end(), s.begin(), s.end());
    sequense.insert(sequense.end(), s12.begin(), s12.end());

    for(int i = 0; i < sequense.size(); i++){
        cout<< sequense[i];
    }
    //vector<char> s = sequense;
    vector<double> p = prob(sequense,alphabet);
    for(int tt = 0; tt < alphabet; tt++){
        cout<<p[tt]<<" ";
    }
    cout<<"\n";
    for(int n = 2; n < 25; n++){
        cout<<"a\n";
        vector<double> vector_of_maxf;
        for(int g = 0; g < 20; g++) {
            vector<vector<vector<int>>> set_of_rand_mat;
            //cout<<"b\n";
            while (set_of_rand_mat.size() < 100) {//даже 10**4 работает хуилиард времени (я ждал 2 минуты), ты че дядя, какие 10**8
                vector<vector<int>> mat = random_matrix(alphabet, n);
                if ((int) set_of_rand_mat.size() == 0) {
                    set_of_rand_mat.push_back(mat);
                } else {
                    double minim = 2;
                    for (const auto &j: set_of_rand_mat) {
                        double k = measure_of_difference(j, mat, alphabet, n);
                        if (k < minim) {
                            minim = k;
                        }
                    }
                    if (minim >= 1.0) {
                        set_of_rand_mat.push_back(mat);
                    }
                }
            }
            cout<<"c\n";
            double maxf = 0.0;
            for (const auto &k: set_of_rand_mat) {
                double B = -1.5;
                vector<vector<double>> mod_mat = random_matrix_modified(k, alphabet, n, p, alphabet * n, B);
                for(int i = 0; i < alphabet; i++){
                    for(int j = 0; j < n; j++){
                        cout<< mod_mat[i][j]<<" | ";
                    }
                    cout<<"\n";

                }
                pair<vector<vector<double>>, vector<vector<int>>> cur_pair = F(mod_mat, 0.6, sequense, n);
                int diff = difference(cur_cout<<"k\n";
                double k2 = 1.5;
                while (diff < 440 or diff > 560){
                    //cout<<"g "<<diff<<"\n";
                    k2 = k2/2;
                    if (diff > 660) {
                        B = B - k2;
                        mod_mat = random_matrix_modified(k, alphabet, n, p, alphabet * n, B);
                        cur_pair = F(mod_mat, 0.6, sequense, n);
                        diff = difference(cur_pair);
                    } else {
                        B = B + k2;
                        mod_mat = random_matrix_modified(k, alphabet, n, p, alphabet * n, B);
                        cur_pair = F(mod_mat, 0.6, sequense, n);
                        diff = difference(cur_pair);
                    }
                }
                cout<<"d\n";
                double f2 = get<2>(maxim(cur_pair.first, (int) cur_pair.first.size()));
                if (maxf < f2) {
                    maxf = f2;
                }
            }
            vector_of_maxf.push_back(maxf);
            sequense = randomize(sequense);
            cout<<"e\n";
        }
        double a = vector_of_maxf[0];
        vector_of_maxf.erase(vector_of_maxf.begin());
        double sigma = findAverage(vector_of_maxf);
        double dispersion = findVariance(vector_of_maxf);
        double element_end = (a - sigma)/ sqrt(dispersion);
        res.push_back(element_end);
        cout<<"f";
    }
    for(double re : res){
        cout<< re <<" ";
    }
    return 0;
}