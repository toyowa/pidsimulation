/*
信号とプロペラ回転数の間にタイムラグがある場合のPID制御
剛体モデル
By Toyoaki Washida 02/25/2022
コンパイルコマンド
g++ -o timelag timelag-2.cpp
*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <time.h>

std::string getDate();

//スロットル（回転数:rpm）	
double V0 = 1.30e4;
//double V0 = 6500;
//回転数を力に変換 後で設定	
double Delta = 1.0e-8;
//中心とモーター間の距離 m	
//double L = 0.215; 
double L = 0.37; 
//制御パラメータP	
double P = 300;
//制御パラメータD	
double D = 300;
//Δt （時間間隔：秒）	
double h = 0.01;
// タイムラグ期間 秒では tau*h
const int tau = 20;
// シミュレーション期間
//const int T = 3000;
const int T = 2000;
// thetaの初期値
double init_theta = 0.5;
// omegaの初期値
double init_omega = 0.0;

// 係数 η
double eta;

// 慣性モーメント
double I = 0.01; 

//保存される系列
double thetas[T+tau];
double omegas[T+tau];

void setParameters(){
    eta = 4*Delta*L*V0/I;
    std::cout << "Parameter eta = " << eta << "\n";
}

double lagparam(double theta, double omega){
    return -4.0*(L*Delta*V0/I)*(P*theta+D*omega);
}

void rkmethod(int t){
    double k_theta_1 = h*omegas[t];
    double k_omega_1 = h*(lagparam(thetas[t-tau], omegas[t-tau]));
    double k_theta_2 = h*(omegas[t]+k_omega_1/2.0);
    double k_omega_2 = h*(lagparam(thetas[t-tau], omegas[t-tau]));
    double k_theta_3 = h*(omegas[t]+k_omega_2/2.0);
    double k_omega_3 = h*(lagparam(thetas[t-tau], omegas[t-tau]));
    double k_theta_4 = h*(omegas[t]+k_omega_3);
    double k_omega_4 = h*(lagparam(thetas[t-tau], omegas[t-tau]));
    thetas[t+1] = thetas[t] + (1.0/6.0)*(k_theta_1+2*k_theta_2+2*k_theta_3+k_theta_4);
    omegas[t+1] = omegas[t] + (1.0/6.0)*(k_omega_1+2*k_omega_2+2*k_omega_3+k_omega_4);
    //std::cout << "lagparam(thetas[t-tau], omegas[t-tau]) = " << lagparam(thetas[t-tau], omegas[t-tau]) << "\n";
    //std::cout << "k_theta_1 = " << k_theta_1 << ", k_theta_2 = " << k_theta_2 
    //          << ", k_theta_3 = " << k_theta_3 << ", k_theta_4 = " << k_theta_4 << "\n";
    //std::cout << "thetas[t] + (1/6)*(k_theta_1+2*k_theta_2+2*k_theta_3+k_theta_4) = " 
    //            << thetas[t] + (1/6)*(k_theta_1+2*k_theta_2+2*k_theta_3+k_theta_4) << "\n";
    //std::cout << "No." << t+1 << "," << thetas[t+1] << "," << omegas[t+1] << "\n";
}

int main(){
    std::cout << "----- 信号とプロペラ回転数の間にタイムラグがある場合のPID制御の実行 -----\n";
    std::cout << "----- 剛体としてみたドローン -----\n";
    setParameters();
    // 初期値をセットする
    for(int i=0;i<=tau;i++){
        thetas[i] = init_theta;
        omegas[i] = init_omega;
        //std::cout << "No." << i << ", theta = " << thetas[i] << ", omega = " << omegas[i] << "\n";
    }
    std::cout << "----- シミュレーションの開始 -----\n";
    for(int i=0;i<T-1;i++){
        rkmethod(i+tau);
    }
    std::cout << "----- シミュレーションの終了 -----\n";

    // 計算結果の書き出し
    std::ofstream file;
    std::string filename = "vflight-" + getDate() + ".txt";
    file.open(filename, std::ios::out);
    file << "----- Set parameters -----\n";
    file << "V0 = " << V0 << "\n";
    file << "Delta = " << Delta << "\n";
    file << "I = " << I << "\n";
    file << "eta = " << eta << "\n";
    file << "L = " << L << "\n";
    file << "P = " << P << "\n";
    file << "D = " << D << "\n";
    file << "h = " << h << "\n";
    file << "T = " << T << "\n";
    file << "tau = " << tau << "\n";
    file << "init_theta = " << init_theta << "\n";
    file << "init_omega = " << init_omega << "\n";
    file << " ---- Results -----\n";
    file << ",No.,Pitcha angle,rpm,angular velocity\n";
    for(int i=1;i<T+tau;i++){
        file << "No.," << i << "," << thetas[i] << "," 
             << V0-P*thetas[i]-D*(thetas[i]-thetas[i-1])/h << ","
             << omegas[i] << ","
             << "\n";
    }
    file.close();
    return 0;
}

std::string getDate(){
    time_t t = time(nullptr);
    const tm* lt = localtime(&t);
    std::stringstream s;
    s<<"20";
    s<<lt->tm_year-100; //100を引くことで20xxのxxの部分になる
    s<<"-";
    s<<lt->tm_mon+1; //月を0からカウントしているため
    s<<"-";
    s<<lt->tm_mday; //そのまま
    s<<"-";
    s<<lt->tm_hour;
    s<<"-";
    s<<lt->tm_min;
    s<<"-";
    s<<lt->tm_sec;
    std::string result = s.str();
    return result;
}
