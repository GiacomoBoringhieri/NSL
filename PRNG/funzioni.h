#ifndef FUNZIONI_H
#define FUNZIONI_H

#include "random.h"

class FunzioneBase {
public:
    virtual double eval(double x) const = 0; // Metodo virtuale puro
    void integraBlocchi(int M, int N, Random& rnd, const std::string& file_name) const; // Metodo condiviso
    virtual ~FunzioneBase() {} // Distruttore virtuale
    void integraBlocchiCDF(int M, int N, Random& rnd, FunzioneBase& P_r, const std::string& file_name) const;
};

class Linear : public FunzioneBase {
public:
    double eval(double x) const override; 
};

class Cos_g : public FunzioneBase {
public:
    double eval(double x) const override; 
};

class Cos : public FunzioneBase {
public:
    double eval(double x) const override;
};

class var : public FunzioneBase{
public:
    double eval(double x) const override;
};

class P_r_2 : public FunzioneBase{
public:
    double eval(double x) const override;
};

class P_r_exp : public FunzioneBase{
public:
    // Costruttore
    P_r_exp(double x) : lambda(x) {}

    double eval(double x) const override;
private:
    double lambda;
};

class P_r_Lorentz : public FunzioneBase{
    public:
        // Costruttore
        P_r_Lorentz(double x0, double gamma) : x0(x0), gamma(gamma) {};

        double eval(double x) const override;
        
    private:
        double x0, gamma;
};

//class Brownian_Motion : public FunzioneBase{
//public:
//        double eval(double x) const override;

//};


#endif // FUNZIONI_H