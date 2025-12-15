#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <string>

using namespace std;

// Моном представлен как вектор степеней переменных
using Monomial = vector<int>;
// Полином представлен как map: моном -> коэффициент
using Polynomial = map<Monomial, double>;

enum Order { LEX, GRLEX, GREVLEX };

// Сравнение мономов
bool compare_monomials(const Monomial& m1, const Monomial& m2, Order order) {
    if (order == LEX) {
        for (size_t i = 0; i < m1.size(); i++) {
            if (m1[i] > m2[i]) return true;
            if (m1[i] < m2[i]) return false;
        }
        return false;
    }

    if (order == GRLEX) {
        int deg1 = accumulate(m1.begin(), m1.end(), 0);
        int deg2 = accumulate(m2.begin(), m2.end(), 0);
        if (deg1 > deg2) return true;
        if (deg1 < deg2) return false;

        for (size_t i = 0; i < m1.size(); i++) {
            if (m1[i] > m2[i]) return true;
            if (m1[i] < m2[i]) return false;
        }
        return false;
    }

    if (order == GREVLEX) {
        int deg1 = accumulate(m1.begin(), m1.end(), 0);
        int deg2 = accumulate(m2.begin(), m2.end(), 0);
        if (deg1 > deg2) return true;
        if (deg1 < deg2) return false;

        for (int i = m1.size() - 1; i >= 0; i--) {
            if (m1[i] < m2[i]) return true;
            if (m1[i] > m2[i]) return false;
        }
        return false;
    }

    return false;
}

// Находит старший член полинома
pair<Monomial, double> leading_term(const Polynomial& poly, Order order) {
    if (poly.empty()) {
        return {{}, 0.0};
    }

    auto it = max_element(poly.begin(), poly.end(),
        [order](const auto& a, const auto& b) {
            return !compare_monomials(a.first, b.first, order);
        });

    return *it;
}

// Умножение монома на моном
Monomial multiply_monomials(const Monomial& m1, const Monomial& m2) {
    Monomial result(m1.size());
    for (size_t i = 0; i < m1.size(); i++) {
        result[i] = m1[i] + m2[i];
    }
    return result;
}

// НОК двух мономов
Monomial lcm_monomials(const Monomial& m1, const Monomial& m2) {
    Monomial result(m1.size());
    for (size_t i = 0; i < m1.size(); i++) {
        result[i] = max(m1[i], m2[i]);
    }
    return result;
}

// Деление монома на моном (если делится)
pair<bool, Monomial> divide_monomials(const Monomial& m1, const Monomial& m2) {
    Monomial result(m1.size());
    for (size_t i = 0; i < m1.size(); i++) {
        if (m1[i] < m2[i]) {
            return {false, {}};
        }
        result[i] = m1[i] - m2[i];
    }
    return {true, result};
}

// Умножение полинома на моном с коэффициентом
Polynomial multiply_poly_monomial(const Polynomial& poly, const Monomial& mon, double coeff) {
    Polynomial result;
    for (const auto& [m, c] : poly) {
        Monomial new_mon = multiply_monomials(m, mon);
        result[new_mon] = c * coeff;
    }
    return result;
}

// Сложение полиномов
Polynomial add_polynomials(const Polynomial& p1, const Polynomial& p2) {
    Polynomial result = p1;
    for (const auto& [m, c] : p2) {
        result[m] += c;
        if (abs(result[m]) < 1e-10) {
            result.erase(m);
        }
    }
    return result;
}

// Вычитание полиномов
Polynomial subtract_polynomials(const Polynomial& p1, const Polynomial& p2) {
    Polynomial result = p1;
    for (const auto& [m, c] : p2) {
        result[m] -= c;
        if (abs(result[m]) < 1e-10) {
            result.erase(m);
        }
    }
    return result;
}

// S-полином
Polynomial s_polynomial(const Polynomial& p, const Polynomial& q, Order order) {
    auto [lm_p, lc_p] = leading_term(p, order);
    auto [lm_q, lc_q] = leading_term(q, order);

    Monomial lcm = lcm_monomials(lm_p, lm_q);

    auto [_, quot_p] = divide_monomials(lcm, lm_p);
    auto [__, quot_q] = divide_monomials(lcm, lm_q);

    Polynomial term1 = multiply_poly_monomial(p, quot_p, 1.0 / lc_p);
    Polynomial term2 = multiply_poly_monomial(q, quot_q, 1.0 / lc_q);

    return subtract_polynomials(term1, term2);
}

// Редукция полинома по базису
Polynomial reduce(Polynomial poly, const vector<Polynomial>& basis, Order order) {
    bool changed = true;

    while (changed && !poly.empty()) {
        changed = false;
        auto [lm_r, lc_r] = leading_term(poly, order);

        for (const auto& g : basis) {
            if (g.empty()) continue;

            auto [lm_g, lc_g] = leading_term(g, order);
            auto [divisible, quotient] = divide_monomials(lm_r, lm_g);

            if (divisible) {
                Polynomial to_subtract = multiply_poly_monomial(g, quotient, lc_r / lc_g);
                poly = subtract_polynomials(poly, to_subtract);
                changed = true;
                break;
            }
        }
    }

    return poly;
}

// Алгоритм Бухбергера
vector<Polynomial> buchberger(const vector<Polynomial>& F, Order order) {
    vector<Polynomial> G;
    for (const auto& f : F) {
        if (!f.empty()) {
            G.push_back(f);
        }
    }

    while (true) {
        vector<Polynomial> G_prev = G;

        for (size_t i = 0; i < G_prev.size(); i++) {
            for (size_t j = i + 1; j < G_prev.size(); j++) {
                Polynomial s = s_polynomial(G_prev[i], G_prev[j], order);
                Polynomial s_red = reduce(s, G_prev, order);

                if (!s_red.empty()) {
                    G.push_back(s_red);
                }
            }
        }

        if (G.size() == G_prev.size()) {
            break;
        }
    }

    return G;
}

// Вывод полинома
void print_polynomial(const Polynomial& poly, const vector<string>& var_names) {
    bool first = true;
    for (auto it = poly.rbegin(); it != poly.rend(); ++it) {
        const auto& [mon, coeff] = *it;

        if (abs(coeff) < 1e-10) continue;

        if (!first && coeff > 0) cout << " + ";
        if (!first && coeff < 0) cout << " - ";
        if (first && coeff < 0) cout << "-";

        first = false;

        double abs_coeff = abs(coeff);
        bool is_one = abs(abs_coeff - 1.0) < 1e-10;
        bool is_constant = all_of(mon.begin(), mon.end(), [](int d) { return d == 0; });

        if (!is_one || is_constant) {
            cout << abs_coeff;
        }

        for (size_t i = 0; i < mon.size(); i++) {
            if (mon[i] > 0) {
                cout << var_names[i];
                if (mon[i] > 1) {
                    cout << "^" << mon[i];
                }
            }
        }
    }
    if (first) cout << "0";
}

int main() {
    // Пример: F = [x^2 + y, xy + 1]
    // Переменные: x (индекс 0), y (индекс 1)
    vector<string> vars = {"x", "y"};

    // x^2 + y: {2,0}->1, {0,1}->1
    Polynomial p1 = {{{2, 0}, 1.0}, {{0, 1}, 1.0}};

    // xy + 1: {1,1}->1, {0,0}->1
    Polynomial p2 = {{{1, 1}, 1.0}, {{0, 0}, 1.0}};

    vector<Polynomial> F = {p1, p2};

    cout << "Входной базис:" << endl;
    for (size_t i = 0; i < F.size(); i++) {
        cout << "f_" << i+1 << " = ";
        print_polynomial(F[i], vars);
        cout << endl;
    }

    cout << "\nLEX:" << endl;
    auto G_lex = buchberger(F, LEX);
    for (size_t i = 0; i < G_lex.size(); i++) {
        cout << "g_" << i+1 << " = ";
        print_polynomial(G_lex[i], vars);
        cout << endl;
    }

    cout << "\nGRLEX:" << endl;
    auto G_grlex = buchberger(F, GRLEX);
    for (size_t i = 0; i < G_grlex.size(); i++) {
        cout << "g_" << i+1 << " = ";
        print_polynomial(G_grlex[i], vars);
        cout << endl;
    }

    cout << "\nGREVLEX:" << endl;
    auto G_grevlex = buchberger(F, GREVLEX);
    for (size_t i = 0; i < G_grevlex.size(); i++) {
        cout << "g_" << i+1 << " = ";
        print_polynomial(G_grevlex[i], vars);
        cout << endl;
    }

    return 0;
}
