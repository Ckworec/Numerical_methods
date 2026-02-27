#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>

size_t f_count = 0;

// ===== Арифметика =====
template<typename T>
T add(const T& a, const T& b) { return a + b; }

template<typename T>
T scale(const T& a, double s) { return a * s; }

std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> r(a.size());
    for (size_t i = 0; i < a.size(); ++i) r[i] = a[i] + b[i];
    return r;
}

std::vector<double> scale(const std::vector<double>& a, double s) {
    std::vector<double> r(a.size());
    for (size_t i = 0; i < a.size(); ++i) r[i] = a[i] * s;
    return r;
}

template<typename T>
double ro(const T& y1, const T& y2, int p) {
    if constexpr (std::is_same_v<T, std::vector<double>>) {
        double m = 0.0;
        for (size_t i = 0; i < y1.size(); ++i) {
            double d = std::abs((y2[i] - y1[i]) / (std::pow(2.0, p) - 1));
            if (d > m) m = d;
        }
        return m;
    } else {
        return (y2 - y1) / (std::pow(2.0, p) - 1);
    }
}

//////////////////////////////////////////////////////////////
// 1️⃣ Метод порядка O(h^2) по фото
//////////////////////////////////////////////////////////////

template<typename Y>
Y first_method(std::function<Y(double,const Y&)> f,
               double a,double b,const Y& y0,
               double eps=1e-7,double h0=0.1,
               int p=2,double K=8.0)
{
    f_count = 0;
    auto fw = [&](double x,const Y& y){ ++f_count; return f(x,y); };

    double h = h0;
    std::vector<double> x = {a};
    std::vector<Y> y = {y0};
    std::vector<Y> yline = {y0};

    size_t n = 0, iter = 0;

    while (x[n] < b - 1e-15) {
        if (x[n] + h > b) h = b - x[n];

        Y k1,k2,k3,y_temp,y_new,yh,yline_new;

        // ----- полный шаг -----
        k1 = scale(fw(x[n],y[n]),h); iter++;

        y_temp = add(y[n], scale(k1,0.5));
        k2 = scale(fw(x[n]+0.5*h,y_temp),h); iter++;

        y_temp = add(add(y[n],scale(k1,-1.0)), scale(k2,2.0));
        k3 = scale(fw(x[n]+h,y_temp),h); iter++;

        y_new = add(y[n],
                    scale(add(add(k1,scale(k2,4.0)),k3),1.0/6.0));

        y.push_back(y_new);

        // ----- два полу-шага -----
        double hh = h/2;

        k1 = scale(fw(x[n],y[n]),hh); iter++;
        y_temp = add(y[n], scale(k1,0.5));
        k2 = scale(fw(x[n]+0.5*hh,y_temp),hh); iter++;
        y_temp = add(add(y[n],scale(k1,-1.0)), scale(k2,2.0));
        k3 = scale(fw(x[n]+hh,y_temp),hh); iter++;
        yh = add(y[n],
                 scale(add(add(k1,scale(k2,4.0)),k3),1.0/6.0));

        k1 = scale(fw(x[n]+hh,yh),hh); iter++;
        y_temp = add(yh, scale(k1,0.5));
        k2 = scale(fw(x[n]+hh+0.5*hh,y_temp),hh); iter++;
        y_temp = add(add(yh,scale(k1,-1.0)), scale(k2,2.0));
        k3 = scale(fw(x[n]+h,y_temp),hh); iter++;
        yline_new = add(yh,
                        scale(add(add(k1,scale(k2,4.0)),k3),1.0/6.0));

        yline.push_back(yline_new);

        double roloc = ro(y[n+1],yline[n+1],p);

        if (std::abs(roloc) > eps) {
            h/=2;
            y.pop_back();
            yline.pop_back();
        } else {
            if constexpr (std::is_same_v<Y,std::vector<double>>) {
                Y corr(y_new.size());
                for(size_t i=0;i<corr.size();++i)
                    corr[i]=yline_new[i]+(yline_new[i]-y_new[i])/(std::pow(2.0,p)-1);
                y[n+1]=corr;
            } else {
                y[n+1]=yline_new+(yline_new-y_new)/(std::pow(2.0,p)-1);
            }

            x.push_back(x[n]+h);
            if(std::abs(roloc)<eps/K)
                h=std::min({h*2,h0*2,b-x[n+1]});
            n++;
        }
    }

    std::cout<<"Iter1="<<iter<<" f="<<f_count<<"\n";
    return y.back();
}

//Метод порядка O(h^3)

template<typename Y>
Y second_method(std::function<Y(double,const Y&)> f,
                double a,double b,const Y& y0,
                double eps=1e-7,double h0=0.1,
                int p=3,double K=16.0)
{
    f_count=0;
    auto fw=[&](double x,const Y& y){++f_count;return f(x,y);};

    double h=h0;
    std::vector<double> x={a};
    std::vector<Y> y={y0};
    std::vector<Y> y2={y0};

    size_t n=0,iter=0;

    while(x[n]<b-1e-15){
        if(x[n]+h>b) h=b-x[n];

        Y k1,k2,k3,k4,y_temp,y_new;

        k1=scale(fw(x[n],y[n]),h);iter++;
        y_temp=add(y[n],scale(k1,0.5));
        k2=scale(fw(x[n]+0.5*h,y_temp),h);iter++;
        y_temp=add(y[n],scale(k2,0.5));
        k3=scale(fw(x[n]+0.5*h,y_temp),h);iter++;
        y_temp=add(y[n],k3);
        k4=scale(fw(x[n]+h,y_temp),h);iter++;

        y_new=add(y[n],
                  scale(add(add(k1,scale(k2,2.0)),
                            add(scale(k3,2.0),k4)),1.0/6.0));

        y.push_back(y_new);
        y2.push_back(y_new);

        x.push_back(x[n]+h);
        n++;
    }

    std::cout<<"Iter2="<<iter<<" f="<<f_count<<"\n";
    return y.back();
}

// Метод с контрольным членом

template<typename Y>
Y third_method(std::function<Y(double,const Y&)> f,
               double a,double b,const Y& y0,
               double eps=1e-4,double h0=0.1,
               int p=3,double K=16.0)
{
    f_count=0;
    auto fw=[&](double x,const Y& y){++f_count;return f(x,y);};

    double h=h0;
    std::vector<double> x={a};
    std::vector<Y> y={y0};

    size_t n=0,iter=0;

    while(x[n]<b-1e-15){
        if(x[n]+h>b) h=b-x[n];

        Y k1,k2,k3,k4,k5,y_temp,y_new;
        double E=0.0;

        k1=scale(fw(x[n],y[n]),h);iter++;

        y_temp=add(y[n],scale(k1,1.0/3.0));
        k2=scale(fw(x[n]+h/3.0,y_temp),h);iter++;

        y_temp=add(y[n],scale(add(k1,k2),1.0/6.0));
        k3=scale(fw(x[n]+h/3.0,y_temp),h);iter++;

        y_temp=add(add(y[n],scale(k1,1.0/8.0)),
                   scale(k3,3.0/8.0));
        k4=scale(fw(x[n]+0.5*h,y_temp),h);iter++;

        y_temp=add(add(add(y[n],scale(k1,0.5)),
                       scale(k3,-1.5)),
                       scale(k4,2.0));
        k5=scale(fw(x[n]+h,y_temp),h);iter++;

        y_new=add(y[n],
                  scale(add(add(k1,scale(k4,4.0)),k5),1.0/6.0));

        if constexpr(std::is_same_v<Y,std::vector<double>>){
            for(size_t i=0;i<y[n].size();++i){
                double Ei=(2.0*k1[i]-9.0*k3[i]
                           +8.0*k4[i]-k5[i])/30.0;
                if(std::abs(Ei)>E) E=std::abs(Ei);
            }
        } else {
            E=std::abs((2.0*k1-9.0*k3+8.0*k4-k5)/30.0);
        }

        if(E>eps){
            h/=2;
        }else{
            y.push_back(y_new);
            x.push_back(x[n]+h);
            if(E<eps/K)
                h=std::min({h*2,h0*2,b-x.back()});
            n++;
        }
    }

    std::cout<<"Iter3="<<iter<<" f="<<f_count<<"\n";
    return y.back();
}

double func(double x,const double&){return x;}

int main(){

    std::vector<double> y0={1.0,1.0};
    auto fvec=[](double,const std::vector<double>& y){
        return std::vector<double>{y[0],y[1]};
    };

    auto r1 = first_method<std::vector<double>>(fvec, 0.0, 1.0, y0);
    std::cout<<"First: ["<<r1[0]<<","<<r1[1]<<"]\n";

    auto r2 = second_method<std::vector<double>>(fvec, 0.0, 1.0, y0);
    std::cout<<"Second: ["<<r2[0]<<","<<r2[1]<<"]\n";

    auto r3 = third_method<std::vector<double>>(fvec, 0.0, 1.0, y0);
    std::cout<<"Third: ["<<r3[0]<<","<<r3[1]<<"]\n";

    return 0;
}