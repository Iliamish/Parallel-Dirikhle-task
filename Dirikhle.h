#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <mutex>
#include <tbb/tbb.h>

const double eps_lim = std::numeric_limits<double>::epsilon();

double equalsZero(double num) {
    return (num < eps_lim&& num > -eps_lim) ? 0.0 : num;
}

double m1(double x, double y) {
    return equalsZero(sin(M_PI * y) * sin(M_PI * y));
}

double m2(double x, double y) {
    return equalsZero(sin(M_PI * 2 * y) * sin(M_PI * 2 * y));
}

double m3(double x, double y) {
    return equalsZero(sin(M_PI * x) * sin(M_PI * x));
}

double m4(double x, double y) {
    return equalsZero(sin(M_PI * 2 * x) * sin(M_PI * 2 * x));
}

double function(double x, double y) {
    return abs(x * x - 2 * y);
}

void initBorders(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d) {
    int i, j;
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < m + 1; j++) {
            v[i].push_back(0);
        }
    }

    for (j = 0; j < m + 1; j++) {
        v[0][j] = m1(a, c + (d - c) / m * j);
    }
    for (j = 0; j < m + 1; j++) {
        v[n][j] = m2(b, c + (d - c) / m * j);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][0] = m3(a + i * (b - a) / n, c);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][m] = m4(a + i * (b - a) / n, d);
    }
}

int simple_version(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d, int Nmax, double eps) {
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    int S =0;
    double eps_max = 0;
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    initBorders(v, n, m, a, b, c, d);

    while (!flag) {
        eps_max = 0;
        for (j = 1; j < m; j++)
            for (i = 1; i < n; i++) {
                v_old = v[i][j];
                v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                double p = a + i * (b - a) / n;
                double e = c + j * (d - c) / m;
                double func_t = function(p, e);
                v_new = v_new + func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                if (eps_cur > eps_max) { eps_max = eps_cur; };
                v[i][j] = v_new;
            }
        S = S + 1;
        if ((eps_max < eps) || (S >= Nmax)) { flag = true; }
    }
    return S;
}

int simple_tbb_version(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d, int Nmax, double eps) {
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    int S = 0;
    double eps_max = 0;
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);



    initBorders(v, n, m, a, b, c, d);

    std::mutex mutex;

    auto lambda = [&](const tbb::blocked_range2d<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_tmp = 0, eps_cur;
        for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
            for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                v_old = v[i][j];
                v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                p = a + i * (b - a) / n;
                e = c + j * (d - c) / m;
                func_t = function(p, e);
                v_new = v_new + func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                eps_tmp = max(eps_cur, eps_tmp);
                v[i][j] = v_new;
            }
        }
        mutex.lock();
        if (eps_cur > eps_max) {
            eps_max = eps_cur;
        };
        mutex.unlock();
    };
        
    while (!flag) {
        eps_max = 0;
        tbb::parallel_for(tbb::blocked_range2d<size_t>(1, m, 1, n), lambda);
        S = S + 1;
        if ((eps_max < eps) || (S >= Nmax)) { flag = true; }
    }
    return S;
}

int tbb_version_2(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d, int Nmax, double eps) {
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    int S = 0;
    std::atomic<double> eps_max = 0;
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    initBorders(v, n, m, a, b, c, d);

    auto lambda = [&](const tbb::blocked_range2d<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_cur;
        for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
            for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                v_old = v[i][j];
                v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                p = a + i * (b - a) / n;
                e = c + j * (d - c) / m;
                func_t = function(p, e);
                v_new = v_new + func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                v[i][j] = v_new;
                if (eps_cur > eps_max.load()) { eps_max.store(eps_cur); };
            }
        }
    };

    while (!flag) {
        eps_max = 0;
        tbb::parallel_for(tbb::blocked_range2d<size_t>(1, m, 1, n), lambda);
        S = S + 1;
        if ((eps_max < eps) || (S >= Nmax)) { flag = true; }
    }
    return S;
}

int tbb_version_2_1(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d, int Nmax, double eps) {
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    int S = 0;
    std::atomic<double> eps_max = 0;
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    initBorders(v, n, m, a, b, c, d);

    auto lambda = [&](const tbb::blocked_range2d<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_cur;
        for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
            for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                v_old = v[i][j];
                v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                p = a + i * (b - a) / n;
                e = c + j * (d - c) / m;
                func_t = function(p, e);
                v_new = v_new + func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                v[i][j] = v_new;
                if (eps_cur > eps_max.load(std::memory_order_acquire)) { eps_max.store(eps_cur, std::memory_order_release); };
            }
        }
    };

    while (!flag) {
        eps_max = 0;
        tbb::parallel_for(tbb::blocked_range2d<size_t>(1, m, 1, n), lambda);
        S = S + 1;
        if ((eps_max < eps) || (S >= Nmax)) { flag = true; }
    }
    return S;
}

int tbb_version_3(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d, int Nmax, double eps) {
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    int S = 0;
    std::atomic<double> eps_max = 0;
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    initBorders(v, n, m, a, b, c, d);

    auto lambda = [&](const tbb::blocked_range2d<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_cur, tmp = 0;
        for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
            for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                v_old = v[i][j];
                v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                p = a + i * (b - a) / n;
                e = c + j * (d - c) / m;
                func_t = function(p, e);
                v_new = v_new + func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                v[i][j] = v_new;
                if (eps_cur > tmp) {tmp = eps_cur; };
            }
        }
        if (tmp > eps_max.load()) { eps_max.store(tmp); };
    };

    while (!flag) {
        eps_max = 0;
        tbb::parallel_for(tbb::blocked_range2d<size_t>(1, m, 1, n), lambda);
        S = S + 1;
        if ((eps_max < eps) || (S >= Nmax)) { flag = true; }
    }
    return S;
}


int tbb_version_4(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d, int Nmax, double eps) {
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    int S = 0;
    std::atomic<double> eps_max = 0;
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    initBorders(v, n, m, a, b, c, d);

    auto lambda = [&](const tbb::blocked_range<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_cur, thread_eps_max;
        for (size_t i = r.begin(); i != r.end(); ++i) {
            thread_eps_max = 0;
            for (size_t j = 1; j != n; ++j) {
                v_old = v[i][j];
                v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                p = a + i * (b - a) / n;
                e = c + j * (d - c) / m;
                func_t = function(p, e);
                v_new = v_new + func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                v[i][j] = v_new;   
                if (eps_cur > thread_eps_max) { thread_eps_max = eps_cur; };
            }
            if (thread_eps_max > eps_max.load()) { eps_max.store(thread_eps_max); };
        }

    };

    while (!flag) {
        eps_max = 0;
        tbb::parallel_for(tbb::blocked_range<size_t>(1, m), lambda);
        S = S + 1;
        if ((eps_max < eps) || (S >= Nmax)) { flag = true; }
    }
    return S;
}


int tbb_version_5(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d, int Nmax, double eps) {
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    int S = 0;
    std::atomic<double> eps_max = 0;
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    initBorders(v, n, m, a, b, c, d);

    std::vector<std::vector<double>> v_n;

    auto lambda = [&](const tbb::blocked_range<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_cur, thread_eps_max;
        for (size_t i = r.begin(); i != r.end(); ++i) {
            thread_eps_max = 0;
            for (size_t j = 1; j != n; ++j) {
                v_old = v[i][j];
                v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                p = a + i * (b - a) / n;
                e = c + j * (d - c) / m;
                func_t = function(p, e);
                v_new = v_new + func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                v_n[i][j] = v_new;
                if (eps_cur > thread_eps_max) { thread_eps_max = eps_cur; };
            }
            if (thread_eps_max > eps_max.load()) { eps_max.store(thread_eps_max); };
        }
    };
        
    while (!flag) {
        eps_max = 0;
        v_n = v;
        tbb::parallel_for(tbb::blocked_range<size_t>(1, m), lambda);
        S = S + 1;
        if ((eps_max < eps) || (S >= Nmax)) { flag = true; }
        v = v_n;
    }
    return S;
}

void setNull(std::vector<std::vector<double>>& v) {
    std::cout << "-----------" << std::endl;
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = 0; j < v.size(); ++j) {
            std::cout << v[i][j] << " ";
            v[i][j] = 0;
        }
        std::cout << std::endl;
    }
}

/*only for mxm matrix*/
int tbb_version_6(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d, int Nmax, double eps) {
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    int S = 0;
    std::atomic<double> eps_max = 0;
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    initBorders(v, n, m, a, b, c, d);

    std::vector<double> err_vec(m + 1);

    int nmax;
    int nmin;

    //std::vector<std::vector<double>> null_vec(v.size());

    //for (size_t i = 0; i != v.size(); ++i) {
    //    for (size_t j = 0; j != v.size(); ++j) {
    //        null_vec[i].push_back(0);
    //    }
    //}

    auto lambda1 = [&](const tbb::blocked_range<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_cur, thread_eps_max = 0;
        int k = 0;
        for (size_t i = r.begin(); i != r.end(); ++i) {
            thread_eps_max = 0;
            j = nmax + 1 - i;
            v_old = v[i][j];

            v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
            p = a + i * (b - a) / n;
            e = c + j * (d - c) / m;
            func_t = function(p, e);
            v_new = v_new + func_t;
            v_new = v_new / a2;
            eps_cur = fabs(v_old - v_new);
            v[i][j] = v_new;
            if (eps_cur > err_vec[i]) {
                err_vec[i] = eps_cur;
                //null_vec[i][j] = 1;
            };
        }
    };

    auto lambda2 = [&](const tbb::blocked_range<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_cur, thread_eps_max = 0;
        int k = 0;
        for (size_t i = r.begin(); i != r.end(); ++i) {
            thread_eps_max = 0;
            j = nmin +  nmax - i - 1;
            v_old = v[i][j];
            
            v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
            p = a + i * (b - a) / n;
            e = c + j * (d - c) / m;
            func_t = function(p, e);
            v_new = v_new + func_t;
            v_new = v_new / a2;
            eps_cur = fabs(v_old - v_new);
            v[i][j] = v_new;
            if (eps_cur > err_vec[i]) {
                err_vec[i] = eps_cur;
                //null_vec[i][j] = 1;
            };
        }
    };
    double tmp = 0;
    while (!flag) {
        tmp = 0;
        for (size_t p = 1; p < m; ++p) {
            nmax = p;
            nmin = 1;
            tbb::parallel_for(tbb::blocked_range<size_t>(1, p+1), lambda1);
            //setNull(null_vec);
        }

        for (size_t p = 2; p <= m; ++p) {
            nmax = m;
            nmin = p;
            tbb::parallel_for(tbb::blocked_range<size_t>(p, m), lambda2);
            //setNull(null_vec);
        }

        S = S + 1;

        for (size_t p = 1; p < m; ++p) {
            if (err_vec[p] > tmp) { 
                tmp = err_vec[p]; 
           };
        }
        //eps_max.store(tmp);
        
        if ((tmp < eps) || (S >= Nmax)) { flag = true; }
    }
    return S;
}


/*only for mxm matrix*/
int tbb_version_8(std::vector<std::vector<double>>& v, const int n, const int m, double a, double b, double c, double d, int Nmax, double eps) {
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    int S = 0;
    std::atomic<double> eps_max = 0;
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    initBorders(v, n, m, a, b, c, d);

    std::vector<double> err_vec(m + 1);

    int nmax;
    int nmin;

    std::vector<std::vector<double>> null_vec(v.size());

    for (size_t i = 0; i != v.size(); ++i) {
        for (size_t j = 0; j != v.size(); ++j) {
            null_vec[i].push_back(0);
        }
    }

    int block_size = 2;

    auto lambda1 = [&](const tbb::blocked_range<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_cur, tmp = 0;
        int k = 0;
        for (size_t q = r.begin(); q < r.end(); ++q) {
            for (size_t i = 1 + (q-1) * block_size; i < q * block_size + 1; ++i) {
                for (size_t j = nmax * block_size + 1 - (q) * block_size; j < nmax * block_size + 1 - (q) * block_size + block_size; ++j) {
                    std::cout << i << "--" << j << "--" << q << "--" << r.end() << std::endl;
                    v_old = v[i][j];
                    null_vec[i][j] = 1;
                    v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                    p = a + i * (b - a) / n;
                    e = c + j * (d - c) / m;
                    func_t = function(p, e);
                    v_new = v_new + func_t;
                    v_new = v_new / a2;
                    eps_cur = fabs(v_old - v_new);
                    tmp = max(eps_cur, tmp);
                    v[i][j] = v_new;
                }
            }
        }
        std::cout << "--------------" << std::endl;
        if (tmp > eps_max.load()) { eps_max.store(tmp); };
    };

    auto lambda2 = [&](const tbb::blocked_range<size_t>& r) {
        double v_old, v_new, p, e, func_t, eps_cur, tmp = 0;
        int k = 0;
        /*need to fix*/
        for (size_t q = r.begin(); q < r.end(); ++q) {
            for (size_t i = 1 + (q - 1) * block_size; i != q * block_size + 1; ++i) {
                for (size_t j = nmin * block_size + nmax * block_size - 1 - (q - 1) * block_size; j != nmin * block_size + nmax * block_size - 1 - (q - 1) * block_size + block_size; ++j) {
                    std::cout << i << "--" << j << "--" << q << "--" << r.end() << std::endl;
                    v_old = v[i][j];
                    null_vec[i][j] = 1;
                    v_new = -(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                    p = a + i * (b - a) / n;
                    e = c + j * (d - c) / m;
                    func_t = function(p, e);
                    v_new = v_new + func_t;
                    v_new = v_new / a2;
                    eps_cur = fabs(v_old - v_new);
                    tmp = max(eps_cur, tmp);
                    v[i][j] = v_new;
                }
            }
        }      
    };
    double tmp = 0;
    while (!flag) {
        tmp = 0;
        for (size_t p = 1; p < m / block_size; ++p) {
            nmax = p;
            nmin = 1;
            tbb::parallel_for(tbb::blocked_range<size_t>(1, p + 1), lambda1);
            setNull(null_vec);
        }

        for (size_t p = 2; p <= m / block_size; ++p) {
            nmax = m / block_size;
            nmin = p;
            tbb::parallel_for(tbb::blocked_range<size_t>(p, m), lambda2);
            setNull(null_vec);
        }

        S = S + 1;

        if ((tmp < eps) || (S >= Nmax)) { flag = true; }
    }
    return S;
}