#include <bits/stdc++.h>
#define REP(i, n) for (int i = 0; (i) < (int)(n); ++ (i))
#define REP3(i, m, n) for (int i = (m); (i) < (int)(n); ++ (i))
#define REP_R(i, n) for (int i = int(n) - 1; (i) >= 0; -- (i))
#define REP3R(i, m, n) for (int i = int(n) - 1; (i) >= (int)(m); -- (i))
#define ALL(x) begin(x), end(x)
#define dump(x) cerr << #x " = " << x << endl
using ll = long long;
using namespace std;
template <class T> using reversed_priority_queue = priority_queue<T, vector<T>, greater<T> >;
template <class T, class U> inline void chmax(T & a, U const & b) { a = max<T>(a, b); }
template <class T, class U> inline void chmin(T & a, U const & b) { a = min<T>(a, b); }
template <typename X, typename T> auto vectors(X x, T a) { return vector<T>(x, a); }
template <typename X, typename Y, typename Z, typename... Zs> auto vectors(X x, Y y, Z z, Zs... zs) { auto cont = vectors(y, z, zs...); return vector<decltype(cont)>(x, cont); }
template <typename T> ostream & operator << (ostream & out, vector<T> const & xs) { REP (i, int(xs.size()) - 1) out << xs[i] << ' '; if (not xs.empty()) out << xs.back(); return out; }

class xor_shift_128 {
public:
    typedef uint32_t result_type;
    xor_shift_128(uint32_t seed = 42) {
        set_seed(seed);
    }
    void set_seed(uint32_t seed) {
        a = seed = 1812433253u * (seed ^ (seed >> 30));
        b = seed = 1812433253u * (seed ^ (seed >> 30)) + 1;
        c = seed = 1812433253u * (seed ^ (seed >> 30)) + 2;
        d = seed = 1812433253u * (seed ^ (seed >> 30)) + 3;
    }
    uint32_t operator() () {
        uint32_t t = (a ^ (a << 11));
        a = b; b = c; c = d;
        return d = (d ^ (d >> 19)) ^ (t ^ (t >> 8));
    }
    static constexpr uint32_t max() { return numeric_limits<result_type>::max(); }
    static constexpr uint32_t min() { return numeric_limits<result_type>::min(); }
private:
    uint32_t a, b, c, d;
};

struct photo_t {
    bool is_vertical;
    vector<int> tags;
};

typedef pair<int, int> slide_t;

vector<int> get_tags(slide_t const & slide, vector<photo_t> const & photos) {
    if (slide.second == -1) {
        return photos[slide.first].tags;
    } else {
        vector<int> tags;
        merge(ALL(photos[slide.first].tags), ALL(photos[slide.second].tags), back_inserter(tags));
        tags.erase(unique(ALL(tags)), tags.end());
        return tags;
    }
}

ll get_score_delta(slide_t const & a, slide_t const & b, vector<photo_t> const & photos) {
    auto a_tags = get_tags(a, photos);
    auto b_tags = get_tags(b, photos);
    int x = 0;
    int y = 0;
    int z = 0;
    while (not a_tags.empty() and not b_tags.empty()) {
        if (a_tags.back() > b_tags.back()) {
            ++ x;
            a_tags.pop_back();
        } else if (a_tags.back() < b_tags.back()) {
            ++ z;
            b_tags.pop_back();
        } else {
            ++ y;
            a_tags.pop_back();
            b_tags.pop_back();
        }
    }
    x += a_tags.size();
    z += b_tags.size();
    return min(x, min(y, z));
}

ll compute_score(vector<slide_t> const & slides, vector<photo_t> const & photos) {
    ll acc = 0;
    REP (i, (int)slides.size() - 1) {
        acc += get_score_delta(slides[i], slides[i + 1], photos);
    }
    return acc;
}

template <class Generator>
vector<slide_t> solve(int n, vector<photo_t> const & photos, Generator & gen) {
    vector<slide_t> slides;

    if (getenv("RESUME")) {
        ifstream ifs(getenv("RESUME"));
        int s; ifs >> s;
        while (s --) {
            int i; ifs >> i;
            char c; ifs.get(c);
            if (c == '\n') {
                slides.emplace_back(i, -1);
                assert (not photos[i].is_vertical);
            } else if (c == ' ') {
                int j; ifs >> j;
                slides.emplace_back(i, j);
                assert (photos[i].is_vertical);
                assert (photos[j].is_vertical);
            } else {
                assert (false);
            }
        }

    } else {  // make an initial state
        int vr = -1;
        REP (i, n) {
            if (photos[i].is_vertical) {
                if (vr == -1) {
                    vr = i;
                } else {
                    slides.emplace_back(vr, i);
                    vr = -1;
                }
            } else {
                slides.emplace_back(i, -1);
            }
        }
        assert (vr == -1);
    }

    int s = slides.size();
    ll score = compute_score(slides, photos);

    vector<slide_t> result = slides;
    ll highscore = score;
    cerr << "[*] highscore = " << highscore << endl;

    constexpr int TIME_LIMIT = 1000;  // msec
    chrono::high_resolution_clock::time_point clock_begin = chrono::high_resolution_clock::now();
    double temperature = 1;
    for (unsigned iteration = 0; ; ++ iteration) {
        if (iteration % 128 == 0) {
            chrono::high_resolution_clock::time_point clock_end = chrono::high_resolution_clock::now();
            temperature = 1 - chrono::duration_cast<chrono::milliseconds>(clock_end - clock_begin).count() / TIME_LIMIT;
            if (temperature < 0) {
                cerr << "[*] iteration = " << iteration << ": done" << endl;
                break;
            }
        }

        int i = uniform_int_distribution<int>(0, s - 1)(gen);
        int j = uniform_int_distribution<int>(0, s - 1)(gen);
        if (i > j) swap(i, j);
        if (i == 0) continue;
        if (j == n - 1) continue;
        if (i + 1 == j) continue;
        if (slides[i].second != -1 and slides[j].second != -1) {
            if (bernoulli_distribution(0.5)(gen)) {
                swap(slides[i].first, slides[i].second);
            }
        }
        ll delta = 0;
        delta -= get_score_delta(slides[i - 1], slides[i], photos);
        delta -= get_score_delta(slides[i], slides[i + 1], photos);
        delta -= get_score_delta(slides[j - 1], slides[j], photos);
        delta -= get_score_delta(slides[j], slides[j + 1], photos);
       // if (slides[i].second != -1 and slides[j].second != -1) {
       //     swap(slides[i].first, slides[j].first);
       // } else {
            swap(slides[i], slides[j]);
        //}
        delta += get_score_delta(slides[i - 1], slides[i], photos);
        delta += get_score_delta(slides[i], slides[i + 1], photos);
        delta += get_score_delta(slides[j - 1], slides[j], photos);
        delta += get_score_delta(slides[j], slides[j + 1], photos);

        constexpr double boltzmann = 3;
        if (delta >= 0 or bernoulli_distribution(exp(boltzmann * delta) * temperature)(gen)) {
        // if (delta >= 0) {
            if (delta < 0) {
                cerr << "[*] iteration = " << iteration << ": delta = " << delta << ": p = " << exp(boltzmann * delta) * temperature << endl;
            }
            score += delta;
            if (highscore < score) {
                highscore = score;
                result = slides;
                cerr << "[*] iteration = " << iteration << ": highscore = " << highscore << endl;
            }
        } else {
         //   if (slides[i].second != -1 and slides[j].second != -1) {
         //       swap(slides[i].first, slides[j].first);
         //   } else {
                swap(slides[i], slides[j]);
            //}
        }
    }

    cerr << "[*] highscore = " << highscore << endl;
    return result;
}

int main() {
    // input
    int n; cin >> n;
    unordered_map<string, int> tags;
    vector<photo_t> photos(n);
    REP (i, n) {
        char c; int m; cin >> c >> m;
        photos[i].is_vertical = (c == 'V');
        photos[i].tags.resize(m);
        REP (j, m) {
            string s; cin >> s;
            if (not tags.count(s)) {
                int size = tags.size();
                tags[s] = size;
            }
            photos[i].tags[j] = tags[s];
        }
        sort(ALL(photos[i].tags));
    }

    // solve
    xor_shift_128 gen;
    vector<slide_t> slides = solve(n, photos, gen);

    // output
    vector<bool> used(n);
    cout << slides.size() << endl;
    for (auto [ i, j ] : slides) {
        if (j == -1) {
            cout << i << endl;
        } else {
            cout << i << ' ' << j << endl;
        }
    }
    return 0;
}
