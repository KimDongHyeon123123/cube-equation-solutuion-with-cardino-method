//참고: https://www.youtube.com/watch?v=q14F6fZf5kc

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#define PI 3.141592653589793238462643383279
_C_double_complex make_complex(double real, double imag) {
    _C_double_complex result;
    result._Val[0] = real; // 실수부 설정
    result._Val[1] = imag; // 허수부 설정
    return result;
}

_C_double_complex cmul(_C_double_complex a, _C_double_complex b) {
    _C_double_complex result;

    double a_real = creal(a);
    double a_imag = cimag(a);
    double b_real = creal(b);
    double b_imag = cimag(b);

    double result_real = a_real * b_real - a_imag * b_imag;
    double result_imag = a_real * b_imag + a_imag * b_real;

    result._Val[0] = result_real;
    result._Val[1] = result_imag;

    return result;
}
_C_double_complex add(_C_double_complex a, _C_double_complex b) {
    _C_double_complex result;

    double a_real = creal(a);
    double a_imag = cimag(a);
    double b_real = creal(b);
    double b_imag = cimag(b);

    double result_real = a_real + b_real;
    double result_imag = a_imag + b_imag;

    result._Val[0] = result_real;
    result._Val[1] = result_imag;

    return result;
}
int main()
{
    /*clock_t start, end;
    double cpu_time_used;*/
    
    double co[4]; _C_double_complex roots[3];
    double a, b, c, d, p, q, discrem, u, v,c1,c2,c3;
    for (int i = 3; i > 0; i--) {
        printf("%d차항의 계수 : ", i);
        if (scanf_s("%lf", &co[3 - i]) != 1) {
            printf("잘못된 입력입니다.\n");
            return 1;
        }
    }
    printf("상수항의 계수 : ");
    if (scanf_s("%lf", &co[3]) != 1) {
        printf("잘못된 입력입니다.\n");
        return 1;
    }

    if (co[0] == 0) {
        printf("3차 방정식이 아닙니다. a가 0입니다.\n");
        return 1;
    }
    //start = clock();
    printf("\n입력한 함수\n");
    printf("f(x) = %.1lfx^3 %+0.1lfx^2 %+0.1lfx %+0.1lf\n", co[0], co[1], co[2], co[3]);
    //계산 t^3 + pt+ q = 0, t = x+b/3a
    a = co[0];
    b = co[1];
    c = co[2];
    d = co[3];
    c1 = b / a;
    c2 = c / a;
    c3 = d / a;
    double b3a = -c1 / 3;
    p = c2 - (c1*c1/(3));
    q = c3 - c1*c2/(3) + 2*c1*c1*c1/(27);
    //printf("b3a : %.6lf     p : %.6lf     q : %.6lf\n\n",b3a, p, q);
    //판별식
    discrem = pow(q,2.0) / 4.0 + pow(p,3.0) / 27.0;

    if (discrem > 1e-6) {
        //printf("판별식 > 0\n");
        u = cbrt(-q / 2 + sqrt(discrem));
        v = -p / (3 * u);
        //printf("discrem : %.6lf\n\nu : %.6lf    v : %.6lf\n", discrem, u, v);
        _C_double_complex tu = make_complex(u, 0);
        _C_double_complex tv = make_complex(v, 0);
        roots[0] = make_complex(u +v + b3a, 0);
        _C_double_complex theta = make_complex(0, 2.0 * PI / 3.0);
        _C_double_complex theta2 = make_complex(0, -2.0 * PI / 3.0);

        _C_double_complex tmp1 = cmul(tu, cexp(theta));
        _C_double_complex tmp2 = cmul(tv, cexp(theta));
        _C_double_complex tmp3 = cmul(tu, cexp(theta2));
        _C_double_complex tmp4 = cmul(tv, cexp(theta2));
        roots[1] = add(tmp1, tmp4);
        roots[1]._Val[0] += b3a;
        roots[2] = add(tmp3, tmp2);
        roots[2]._Val[0] += b3a;
    }
    else if (fabs(discrem) < 1e-6) {
        //printf("판별식 == 0\n");
        if (q == 0) { // 3중근
            double root3 = -b / (3 * a);
            roots[0] = make_complex(root3, 0);
            roots[1] = make_complex(root3, 0);
            roots[2] = make_complex(root3, 0);
        }
        else {
            u = cbrt(-q / 2);
            printf("u = %.6lf\n", u);
            roots[0] = make_complex(2.0 * u + b3a, 0);
            roots[1] = make_complex(-u + b3a, 0);
            roots[2] = make_complex(-u + b3a, 0);
        }
    }
    else {//서로다른 3실근
        //printf("판별식 < 0\n");
        //q == 0이면 t(t^2 +p)이므로 x = -b/3a or x = -b/3a +- sqrt(p)
        //q == 0이면 p < 0 왜냐하면 discrem < 0 이기 때문에
        if (q == 0) {
            roots[0] = make_complex(b3a, 0);
            roots[1] = make_complex(b3a + sqrt(-p), 0);
            roots[2] = make_complex(b3a - sqrt(-p), 0);
        }
        else {
            double imgpart = sqrt(-discrem);
            double realpart = -q / 2;
            double r = cbrt(sqrt(imgpart * imgpart + realpart * realpart));//크기
            //printf("r : %.6lf\n     ", r);
            //u + v = t, x = t - b/3a
            for (int n = 0; n < 3; n++)
                roots[n] = make_complex(2 * r * cos(PI / 12.0 + 2.0 * PI * n / 3.0) + b3a, 0);
        }
    }
    printf("\n값\n");
    if (discrem > 1e-6) {
        printf("x1 = %.2lf\n", roots[0]._Val[0]);
        printf("x2 = %.2lf %+0.2lfi\n", roots[1]._Val[0], roots[1]._Val[1]);
        printf("x3 = %.2lf %+0.2lfi\n", roots[2]._Val[0], roots[2]._Val[1]);
    }
    else {
        for (int i = 0; i < 3; i++) {
            printf("x%d= %.2lf\n", i+1, roots[i]._Val[0]);
        }
    }

    //// 실행 시간 측정 종료
    //end = clock();

    //// 시간 계산
    //cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    //printf("실행 시간: %f 초\n", cpu_time_used);
}
