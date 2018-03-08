void getQuatricRoots(double a, double, b, double c, double, d)
{
    double p = (8*a*c - 3*b*b)/(8*a*a);
    double q = (b*b*b -4*a*b*c+8*a*a*d)/(8*a*a*a);

    double R0 = c*c -3*b*d + 12*a*e;
    double R1 = 2*c*c*c -9 *b*c*d +27*b*b*e + 27*a*d*d -72*a*c*e;
    double Q = power((R1+sqrt(R1*R1 - 4*R0*R0*R0))/2, 1/3);
    double S = 0.5*(sqrt((-2/3)*p + 1/(3*a)*(Q + R0/Q)));
    double S2 = S*S;
    double Z = -4*S2 -2*p;
    double Y = -b/(4*a);
    double root;
    if( (Z + q/S) > 0)
    {
        double Imag1 =  0.5 * sqrt(Z +q/S);
        if( Y -S + Imag1>0)
            root = Y -S + Imag1;
        else
            root = Y -S - Imag1;
    }
    else
    {
        double Imag2 =  0.5 * sqrt(Z -q/S);
        if(Y +S + Imag2 > 0)
            root = Y +S + Imag2;
        else
            root = Y +S - Imag2;
    }


}
