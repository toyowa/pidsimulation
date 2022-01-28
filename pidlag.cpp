/*
信号とプロペラ回転数の間にタイムラグがある場合のPID制御
By Toyoaki Washida 01/27/2022
コンパイルコマンド
g++ -o pidlag pidlag.cpp
*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <time.h>

std::string getDate();

double a,b,c,d,e,j,s; // ループでjは使わないこと
//スロットル（回転数:rpm）	
double V0 = 1.30e4;
//回転数を力に変換 後で設定	
double Delta;
//機体質量 Kg	
double m =  0.45;
//重力加速度 m/s~2	
double g = 9.8;
//中心とモーター間の距離 m	
double L = 0.215; 
//制御パラメータP	
double P = 10000;
//制御パラメータD	
double D = 100;
//Δt （時間間隔：秒）	
double h = 0.01;
// タイムラグ期間 秒では tau*h
const int tau = 1;
// シミュレーション期間
//const int T = 3000;
const int T = 1000;
// thetaの初期値
double init_theta = 0.5;
// omegaの初期値
double init_omega = 0.0;

//保存される系列
double thetas[T+tau];
double omegas[T+tau];

void setParameters(){
    // 平衡状態を再現するパラメータ
    std::cout << "基礎パラメータをセットします\n";
    Delta = m*g/(V0*V0);
    //
    a = Delta*V0*V0/(m*L);
    //std::cout << "a = " << a << "\n";
    b = -2*Delta*V0*P/(m*L);
    //std::cout << "b = " << b << "\n";
    c = Delta*P*P/(m*L);
    //std::cout << "c = " << c << "\n";
    d = -2*D*Delta*V0/(m*L);
    //std::cout << "d = " << d << "\n";
    e = 2*D*P*Delta/(m*L);
    //std::cout << "e = " << e << "\n";
    j = Delta*D*D/(m*L);
    //std::cout << "j = " << j << "\n";
    s = -g/L;
    //std::cout << "s = " << s << "\n";
}

double lagparam(double theta, double omega){
    return a+b*theta+c*theta*theta+d*omega+e*theta*omega+j*omega*omega;
}

void rkmethod(int t){
    double k_theta_1 = h*omegas[t];
    double k_omega_1 = h*(lagparam(thetas[t-tau], omegas[t-tau]) + s*cos(thetas[t]));
    double k_theta_2 = h*(omegas[t]+k_omega_1/2.0);
    double k_omega_2 = h*(lagparam(thetas[t-tau], omegas[t-tau]) + s*cos(thetas[t]+k_theta_1/2.0));
    double k_theta_3 = h*(omegas[t]+k_omega_2/2.0);
    double k_omega_3 = h*(lagparam(thetas[t-tau], omegas[t-tau]) + s*cos(thetas[t]+k_theta_2/2.0));
    double k_theta_4 = h*(omegas[t]+k_omega_3);
    double k_omega_4 = h*(lagparam(thetas[t-tau], omegas[t-tau]) + s*cos(thetas[t]+k_theta_3));
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
    file << "m = " << m << "\n";
    file << "g = " << g << "\n";
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
