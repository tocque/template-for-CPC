#include "header.h"

ll exgcd(ll a, ll b, ll &x, ll &y) {
    if(!b) {
        x = 1, y = 0;
        return a;
    }
    ll d = exgcd(b, a%b, y, x);
    y -= a/b*x;
    return d;
}


ll qpow(ll a, ll b, ll p) {
	ll res = 1;
	while(b) {
		if(b&1) res = res*a % p;
		a = a*a % p, b >>= 1;
	}
	return res;
}

// 逆元
struct Inv {
    // 递推求逆元
	const static int N = 1e5;
	ll dat[N];
	void init(int p) {
		dat[1] = 1;
		repr(i, 2, p) dat[i] = (p - p / i) * dat[p % i] % p;
	}
	ll OP([], int v) { return dat[v]; }
    // 欧拉求逆元
    ll OP((), ll a, ll m) {
        static ll x,y;
        if(exgcd(a,m,x,y) == 1) return (x%m + m)%m;
        return -1;
    }
} inv;

/**
 * 
 * 
*/
namespace exlucas {
    ll C(ll n, ll m, ll p) {
        if(m > n) return 0;
        ll res = 1;
        repl(i, 1, m+1, a) {
            a = (n-i) % p;
            res = res*a%p*inv((i+1)%p, p)%p;
        }
        return res;
    }

    ll lucas(ll n, ll m, ll p) {
        return m ? lucas(n/p, m/p, p) * C(n%p, m%p, p) % p : 1;
    }

    ll cal(ll n, ll a, ll b, ll p) {
        if(!n) return 1;
        ll y = n/p, tmp = 1;
        repl(i, 1, p+1) if(i%a) tmp = tmp*i%p;
        ll ans = qpow(tmp, y, p);
        repl(i, y*p+1, n+1) if(i%a) ans = ans*i%p;
        return ans * cal(n/a, a, b, p)%p;
    }

    ll multiLucas(ll n, ll m, ll a, ll b, ll p) {
        ll i, t1, t2, t3, s = 0, tmp;
        for(i = n; i; i/=a) s += i/a;
        for(i = m; i; i/=a) s -= i/a;
        for(i = n-m; i; i/=a) s -= i/a;
        tmp = qpow(a, s, p);
        t1 = cal(n, a, b, p);
        t2 = cal(m, a, b, p);
        t3 = cal(n-m, a, b, p);
        return tmp*t1%p*inv(t2, p)%p*inv(t3, p)%p;
    }


    ll exLucas(ll n, ll m, ll p) {
        static ll q[100], a[100];
        ll t, e = 0;
        for(ll i = 2; i*i <= p; i++) {
            if(!(p%i)) {
                q[++e] = 1;
                t = 0;
                while(!(p%i)) p /= i, q[e] *= i, t++;
                if(q[e] == i) a[e] = lucas(n, m, q[e]);
                else a[e] = multiLucas(n, m, i, t, q[e]);
            }
        }
        if(p > 1) {
            q[++e] = p;
            a[e] = lucas(n, m, p);
        }
        repl(i, 2, e+1, d, c, x, y) {
            d = exgcd(q[1], q[i], x, y);
            c = a[i]-a[1];
            if(c%d) exit(-1);
            t = q[i]/d;
            x = (c/d*x%t+t)%t;
            a[1] = q[1]*x+a[1];
            q[1] = q[1]*q[i]/d;
        }
        return a[1];
    }
}

namespace mobius {
    const ll MOD = 1e9+7, N = 1100000;
    bool tag[N];
    ll cnt, p[N], mob[N];
    ll f[N], g[N], f1[N];
    void init() {
        cnt = 0;
        mob[1] = 1;
        repl(i, 2, N) {
            if(!tag[i]) {
                p[cnt++] = i;
                mob[i] = -1;
            }
            for(ll j = 0; j < cnt && p[j] * i < N; j++) {
                tag[i*p[j]] = 1;
                if(!(i % p[j])) {
                    mob[i*p[j]] = 0;
                    break;
                }
                mob[i*p[j]] = -mob[i];
            }
        }
        f[0] = 0, f[1] = 1;
        repr(i, 2, N) f[i] = (f[i-1] + f[i-2]) % MOD;
        rep(i, N) f1[i] = qpow(f[i], MOD-2, MOD);
        for(ll i = 1; i <= 1000000; i++) g[i] = 1;
        for(ll i = 1; i <= 1000000; i++) {
            for(ll j = i; j <= 1000000; j+=i) {
                ll tmp = mob[j/i];
                if(tmp == -1)
                    umul(g[j], f1[i]);
                else if(tmp == 1)
                    umul(g[j], f[i]);
            }
        }
        g[0] = 1;
        for(ll i = 1; i <= 1000000; i++)
            g[i] = g[i-1] * g[i], g[i] %= MOD;
    }
    void umul(ll &a, ll b) {
        a *= b;
        a %= MOD;
    }
    ll tmp[N];
    void solve(ll n, ll m) {
        ll q = min(n, m);
        ll ans = 1;
        fulset(tmp, 0);
        for(ll i = 1, j = 0; i <= q; i=j+1) {
            ll v = min(n/(n/i), m/(m/i));
            ll pp = (n/i)*(m/i);
            pp %= (MOD-1);
            ans = (ans * qpow(g[v]*inv(g[j], MOD)%MOD, pp, MOD)) % MOD;
            j = v;
        }
        printf("%lld\n", ans);
    }
};

// ------------ 原根
ll root(ll m) {
	static ll v[100000];
    ll buf = m - 1, cnt = 0;
	rep(i, sqrt(buf)+1) {
        if(!(buf % i)) {
            v[cnt++] = i;
            while(!(buf % i)) buf /= i;
        }
    }
    if(buf != 1) v[cnt++] = buf;
    for(int i = 2; ; i++) {
        for(int j = 0; qpow(i, (m-1)/v[j], m) != 1; j++)
        	if(j == cnt) return i;
    }
}

ll euler(ll n) {
    ll ans = n;
    for(int i = 2; i*i<=n; i++) {
        if(!(n%i)) {
            ans = ans / i * (i-1);
            while(!(n%i)) n /= i;
        }
    }
    return (n>1 ? (ans/n*(n-1)) : ans);
}


// ---------------FFT

namespace FFT {
	const int N = 270000;
	const double PI = acos(-1);
	struct cp {
		double x,y;
		cp OP (+, cp&a) {return {x+a.x, y+a.y};}
		cp OP (-, cp&a) {return {x-a.x, y-a.y};}
		cp OP (*, cp&a) {return {x*a.x-y*a.y, x*a.y+y*a.x};}
	} m[N], w[N];
	int r[N];
	void DFT(cp *a, int n, int wp) {
	    rep(i,n) if(i<r[i]) swap(a[i], a[r[i]]);
	    for(int i = 1; i < n; i<<=1) {
	        cp wn = {cos(PI/i), wp*sin(PI/i)};
	        for(int j = i-2; j >= 0; j -= 2) w[j+1] = wn*(w[j] = w[j>>1]);    
			for(int j = 0; j < n; j += i<<1) {
				cp *p = a+j, *q = p+i, x;
				rep(k, i) {
					x = w[k]*q[k];
					q[k] = p[k]-x, p[k] = p[k]+x;
				}
			}
	    }
    }
	void solve (int *a, int an, int *b, int bn, int *c) {
		int n = 1 << int(log2(an+bn)+1);
	    w[0] = {1,0};
		rep(i,en,en = n>>1) r[i<<1|1] = (r[i<<1] = r[i]>>1) + en;
		rep(i,an) m[i].x = a[i];
		rep(i,bn) m[i].y = b[i];
		DFT(m,n,1);
		rep(i,n) m[i] = m[i]*m[i];
		DFT(m,n,-1);
		rep(i,n) c[i] = m[i].y/n/2+0.5;
	}
}

// ------------------ 矩阵
struct matrix {
    const static int N = 305;
	int m[N][N];
	void init() {    //初始化为单位矩阵
		fulset(m, 0);
		rep(i, N) m[i][i] = 1;
	}
	matrix OP (+, matrix& x) {
		matrix r;
		rep(i, N) {
			rep(j, N) r.m[i][j] = (m[i][j] + x.m[i][j]) % MOD;
		}
		return r; 
	}
	matrix OP (*, matrix& x) {
		matrix r;
		rep(i, N) {
			rep(j, N) {
				r.m[i][j] = 0;
				rep(k, N) {
					r.m[i][j] += m[i][k] * x.m[k][j];
					r.m[i][j] %= MOD;
				}
			}
		}
		return r;
	}
	//矩阵快速幂 
	matrix OP (^, int n) {
		matrix r, a = *this;
		r.init();
		while(n) {
			if(n&1) r = r * a;
			a = a * a;
			n >>= 1;
		}
		return r;
	}
}; //矩阵乘法
