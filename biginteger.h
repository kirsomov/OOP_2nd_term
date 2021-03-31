
#include <iostream>
#include<vector>
#include<string>

class BigInteger {
public:
    BigInteger();
    BigInteger(const int& number);
    BigInteger(const std::string& string);

    const std::string toString() const;

    bool Sign() const;
    explicit operator bool() const;

    friend std::istream& operator>> (std::istream& in, BigInteger& number);

    BigInteger& DischargeShift (int shift);

    BigInteger operator- () const;

    BigInteger& operator+= (const BigInteger& number);
    BigInteger& operator-= (const BigInteger& number);
    BigInteger& operator*= (const BigInteger& number);
    BigInteger& operator/= (const BigInteger& number);
    BigInteger& operator%= (const BigInteger& number);

    BigInteger& operator++ ();
    BigInteger& operator-- ();
    const BigInteger operator++ (int);
    const BigInteger operator-- (int);
private:
    bool sign = true; // 0 отрицательный, 1 неотрицательный
    std::vector<unsigned int> digits;
    unsigned int base = 1000000000; //1e9

    unsigned int GetDigitCount() const;
    unsigned int GetDigit(const unsigned int& index) const;

    void RemoveLeadingZeros() ;
    static std::string UIntToString(unsigned int num, const int& min_len);

    friend bool ComparingOfModules(const BigInteger& first, const BigInteger& second);
    BigInteger& AddToModule(const BigInteger& number);
    // корректно работает при вычитании из большего по модулю числа меньшего
    const BigInteger SubtractionFromModule(const BigInteger& number);

    BigInteger& MulToShort(unsigned int number) {
        unsigned int carry = 0;
        for (size_t i = 0; i < digits.size(); i++) {
            long long mul = number * 1ll * digits[i];
            mul += carry;
            digits[i] = mul % base;
            carry = mul / base;
        }
        if (carry > 0)
            digits.push_back(carry);
        RemoveLeadingZeros();
        return *this;
    }

    struct PartsOfNumber;
};

std::ostream& operator<< (std::ostream& out, const BigInteger& number);


bool operator< (const BigInteger& first, const BigInteger& second);
bool operator> (const BigInteger& first, const BigInteger& second);
bool operator<= (const BigInteger& first, const BigInteger& second);
bool operator>= (const BigInteger& first, const BigInteger& second);
bool operator!= (const BigInteger& first, const BigInteger& second);
bool operator== (const BigInteger& first, const BigInteger& second);

const BigInteger operator+ (const BigInteger& first, const BigInteger& second);
const BigInteger operator- (const BigInteger& first, const BigInteger& second);
const BigInteger operator* (const BigInteger& first, const BigInteger& second);
const BigInteger operator/ (const BigInteger& first, const BigInteger& second);
const BigInteger operator% (const BigInteger& first, const BigInteger& second);

class Rational {
public:
    Rational() = default;
    Rational(const int& number);
    Rational(const BigInteger& number);

    const std::string toString() const;
    const std::string asDecimal(size_t precision) const;

    explicit operator double() const;

    Rational operator- () const;

    friend bool operator< (const Rational& first, const Rational& second);

    Rational& operator+= (const Rational& number);
    Rational& operator-= (const Rational& number);
    Rational& operator*= (const Rational& number);
    Rational& operator/= (const Rational& number);


private:
    BigInteger numerator = 0;
    BigInteger denominator = 1;

    const BigInteger GCD(BigInteger first, BigInteger second) const;

    void ReductionToIrreducibleForm();

};

bool operator> (const Rational& first, const Rational& second);
bool operator<= (const Rational& first, const Rational& second);
bool operator>= (const Rational& first, const Rational& second);
bool operator!= (const Rational& first, const Rational& second);
bool operator== (const Rational& first, const Rational& second);

const Rational operator+ (const Rational& first, const Rational& second);
const Rational operator- (const Rational& first, const Rational& second);
const Rational operator* (const Rational& first, const Rational& second);
const Rational operator/ (const Rational& first, const Rational& second);

BigInteger::BigInteger() {
    digits.resize(1);
}

BigInteger::BigInteger(const int& number) {
    sign = number >= 0;
    digits.push_back(abs(number)%base);
    if (abs(number) / base > 0)
        digits.push_back(abs(number) / base);
}

//
BigInteger::BigInteger(const std::string& str) {
    if (str.empty()) {
        sign = true;
        digits.resize(1);
        return;
    }
    sign = str[0] != '-';
    for (int i = static_cast<int>(str.size())-1; i >= 1-sign; i -= 9) {
        unsigned int digit = 0;
        unsigned int ten_pow = 1;
        for (int j = i; j >= std::max(0, i-8); --j) {
            if (j == 0 && !sign)
                break;
            digit += (str[j] - '0') * ten_pow;
            ten_pow *= 10;
        }
        digits.push_back(digit);
    }
}

unsigned int BigInteger::GetDigitCount() const {
    return digits.size();
}

unsigned int BigInteger::GetDigit(const unsigned int& index) const {
    return index < GetDigitCount() ? digits[index] : 0;
}


bool BigInteger::Sign() const{
    return sign;
}

BigInteger::operator bool() const {
    return *this != 0;
}


std::string BigInteger::UIntToString(unsigned int num, const int& min_len) {
    std::string str;
    while (num > 0) {
        str += num % 10 + '0';
        num /= 10;
    }
    int sz = str.size();
    for (int i = 0; i < min_len - sz; i++)
        str += '0';
    sz = str.size();
    for (int i = 0; i < sz / 2; i++)
        std::swap(str[i], str[sz-1-i]);
    return str;
}

void BigInteger::RemoveLeadingZeros()  {
    while(digits[digits.size() - 1] == 0 && digits.size() > 1)
        digits.pop_back();
}

const std::string BigInteger::toString() const {
    int digit_count = GetDigitCount();
    if (digit_count == 0) {
        return "0";
    }
    std::string str;
    if (!sign)
        str += '-';
    str += UIntToString(GetDigit(digit_count-1), 1);
    for (int i = digit_count - 2; i >= 0; i--) {
        unsigned int digit = GetDigit(i);
        str += UIntToString(digit, 9);
    }
    return str;
}

std::ostream& operator<< (std::ostream& out, const BigInteger& number) {
    out << number.toString();
    return out;
}

std::istream& operator>> (std::istream& in, BigInteger& number) {
    std::string str;
    in >> str;
    number = BigInteger(str);
    return in;
}

bool ComparingOfModules(const BigInteger& first, const BigInteger& second) {
    if (first.GetDigitCount() < second.GetDigitCount())
        return true;
    for (int i = first.GetDigitCount() - 1; i >= 0; --i) {
        if (first.GetDigit(i) < second.GetDigit(i))
            return true;
        else
        if (second.GetDigit(i) < first.GetDigit(i))
            return false;
    }
    return false;
}

bool operator< (const BigInteger& first, const BigInteger& second) {
    if (first.Sign() < second.Sign())
        return true;
    if (first.Sign() > second.Sign())
        return false;
    if (first.Sign() && second.Sign())
        return ComparingOfModules(first, second);
    else
        return ComparingOfModules(second, first);
}

bool operator> (const BigInteger& first, const BigInteger& second) {
    return second < first;
}

bool operator<= (const BigInteger& first, const BigInteger& second) {
    return !(second < first);
}

bool operator>= (const BigInteger& first, const BigInteger& second) {
    return !(first < second);
}

bool operator!= (const BigInteger& first, const BigInteger& second) {
    return first < second || second < first;
}

bool operator== (const BigInteger& first, const BigInteger& second) {
    return !(first != second);
}


BigInteger BigInteger::operator- () const {
    BigInteger number(*this);
    number.sign = !number.sign || *this == 0;
    return number;
}

BigInteger& BigInteger::AddToModule(const BigInteger& number) {
    unsigned int carry = 0;
    if (this->GetDigitCount() < number.GetDigitCount())
        this->digits.resize(number.GetDigitCount());
    for (size_t i = 0; i < this->GetDigitCount(); i++) {
        unsigned int sum_of_digits = this->GetDigit(i) + number.GetDigit(i) + carry;
        carry = sum_of_digits / base;
        this->digits[i] = sum_of_digits % base;
        if (i >= number.GetDigitCount() && carry == 0)
            break;
    }
    if (carry > 0)
        this->digits.push_back(carry);
    return *this;
}

const BigInteger BigInteger::SubtractionFromModule(const BigInteger& number) {
    unsigned int carry = 0;
    for (unsigned int i = 0; i < this->GetDigitCount(); i++) {
        int digit = this->GetDigit(i);
        if (i < number.GetDigitCount())
            digit -= number.GetDigit(i);
        digit -= carry;
        carry = 0;
        if (digit < 0) {
            digit += base;
            carry++;
        }
        this->digits[i] = digit;
        if (i >= number.GetDigitCount() && carry == 0)
            break;
    }
    this->RemoveLeadingZeros();
    if (this->digits.size() == 1 && this->digits[0] == 0)
        this->sign = true;
    return *this;
}

BigInteger& BigInteger::operator+= (const BigInteger& number) {
    if (this->Sign() == number.Sign())
        return this->AddToModule(number);
    else
    if (ComparingOfModules(*this, number)) {
        BigInteger temp = *this;
        *this = number;
        this->SubtractionFromModule(temp);
    } else
        this->SubtractionFromModule(number);
    if (this->GetDigit(0) == 0 && this->GetDigitCount() == 1)
        this->sign = true;
    return *this;
}

BigInteger& BigInteger::operator-= (const BigInteger& number) {
    if (this->Sign() != number.Sign())
        return this->AddToModule(number);
    else
    if (ComparingOfModules(*this, number)) {
        BigInteger temp = number;
        *this = temp.SubtractionFromModule(*this);
        this->sign = !this->sign;
    } else
        this->SubtractionFromModule(number);
    if (this->GetDigit(0) == 0 && this->GetDigitCount() == 1)
        this->sign = true;
    return *this;
}

struct BigInteger::PartsOfNumber {
    BigInteger first;
    BigInteger second;

    PartsOfNumber(const BigInteger& number, const unsigned int& part_size) {
        if (part_size < number.GetDigitCount()) {
            std::string str = (number.Sign() ? number : -number).toString();
            std::string first_part = str.substr(0, str.size() - 9 * part_size);
            std::string second_part = str.substr(str.size() - 9 * part_size);
            first = BigInteger(first_part);
            second = BigInteger(second_part);
            first.RemoveLeadingZeros();
            second.RemoveLeadingZeros();
        } else {
            second = (number.Sign() ? number : -number);
            first = 0;
        }
    }
};

BigInteger& BigInteger::DischargeShift (int shift) {
    if (shift >= 0) {
        std::vector<unsigned int> new_digits(shift);
        for (unsigned int digit : digits)
            new_digits.push_back(digit);
        digits = new_digits;
    } else {
        std::vector<unsigned int> new_digits;
        for (int i = -shift; i < static_cast<int>(digits.size()); i++)
            new_digits.push_back(digits[i]);
        if (new_digits.size() == 0)
            new_digits.push_back(0);
        digits = new_digits;
        RemoveLeadingZeros();
    }
    return *this;
}

BigInteger& BigInteger::operator*= (const BigInteger& number) {
    if (this->GetDigitCount() == 1 && number.GetDigitCount() == 1) {
        BigInteger mul;
        unsigned long long digits_mul = this->GetDigit(0);
        digits_mul *= number.GetDigit(0);
        mul.digits[0] = digits_mul%base;
        mul.digits.push_back(digits_mul/base);
        mul.RemoveLeadingZeros();
        mul.sign = this->Sign() == number.Sign() || mul == 0;
        *this = mul;
        return *this;
    }
    if (number.GetDigitCount() == 1 && number.Sign()) {
        this->MulToShort(number.GetDigit(0));
        return *this;
    }
    unsigned int part_size = 0;
    while (true) {
        if (part_size * 2 + 2 <= std::max(this->GetDigitCount(), number.GetDigitCount()))
            part_size++;
        else
            break;
    }
    PartsOfNumber first(*this, part_size);
    PartsOfNumber second(number, part_size);
    // A == first.first
    // B == first.second
    // C == second.first
    // D == second.second
    // (Ax+B)*(Cx+D) = ACx^2+(AD+BC)x+BD
    // first_mul == AC
    // second_mul == BD
    // third_mul == (A+B)(C+D)
    // third_mul-AC-BD == AD+BC
    BigInteger first_mul = first.first;
    first_mul *= second.first;
    BigInteger second_mul = first.second;
    second_mul *= second.second;
    BigInteger third_mul = first.first + first.second;
    third_mul *= second.first + second.second;
    third_mul -= first_mul;
    third_mul -= second_mul;
    first_mul.DischargeShift(part_size*2);
    third_mul.DischargeShift(part_size);
    first_mul.RemoveLeadingZeros();
    second_mul.RemoveLeadingZeros();
    third_mul.RemoveLeadingZeros();
    BigInteger mul = first_mul + second_mul + third_mul;
    mul.RemoveLeadingZeros();
    mul.sign = this->Sign() == number.Sign() || mul == 0;
    *this = mul;
    return *this;
}
/*
BigInteger& BigInteger::operator/= (const BigInteger& number) {
    BigInteger quotient = 0;
    quotient.DischargeShift(this->GetDigitCount()-1);
    for (int i = this->GetDigitCount() - 1; i >= 0; i--) {
        BigInteger current_digit;
        unsigned int left = 0, right = base - 1;
        unsigned int&  mid = quotient.digits[i];
        while (left <= right) {
            mid = (left + right) / 2;
            if (ComparingOfModules(*this, quotient * number))
                right = mid - 1;
            else
                left = mid + 1;
        }
        mid = right;
    }
    quotient.RemoveLeadingZeros();
    quotient.sign = (this->Sign() == number.Sign()) || quotient == 0;
    *this = quotient;
    return *this;
}
*/



BigInteger& BigInteger::operator/= (const BigInteger& number) {
    BigInteger copy = *this;
    copy.sign = true;
    BigInteger module_of_number = number.Sign() ? number : -number;
    for (int i = copy.GetDigitCount() - 1; i >= 0; i--) {
        BigInteger current_factor = 1;
        current_factor.DischargeShift(i);
        unsigned int left = 0, right = base-1;
        BigInteger shifted_copy = copy;
        shifted_copy.DischargeShift(-i);
        unsigned int mid;
        while (left <= right) {
            mid = (left + right) / 2;
            if (shifted_copy < module_of_number * mid)
                right = mid - 1;
            else
                left = mid + 1;
        }
        mid = right;
        this->digits[i] = mid;
        BigInteger mul = module_of_number;
        mul.MulToShort(mid);
        copy -= mul.DischargeShift(i);
    }
    this->RemoveLeadingZeros();
    this->sign = (this->sign == number.Sign());
    if (this->GetDigit(0) == 0 && this->GetDigitCount() == 1)
        this->sign = true;
    return *this;
}

BigInteger& BigInteger::operator%= (const BigInteger& number)  {
    *this -= (*this / number) * number;
    this->RemoveLeadingZeros();
    if (this->GetDigit(0) == 0 && this->GetDigitCount() == 1)
        this->sign = true;
    return *this;
}

const BigInteger operator+ (const BigInteger& first, const BigInteger& second) {
    BigInteger sum = first;
    sum += second;
    return sum;
}

const BigInteger operator- (const BigInteger& first, const BigInteger& second) {
    BigInteger diff = first;
    diff -= second;
    return diff;
}

const BigInteger operator* (const BigInteger& first, const BigInteger& second) {
    BigInteger mul = first;
    mul *= second;
    return mul;
}

const BigInteger operator/ (const BigInteger& first, const BigInteger& second) {
    BigInteger quotient = first;
    quotient /= second;
    return quotient;
}

const BigInteger operator% (const BigInteger& first, const BigInteger& second) {
    BigInteger remainder = first;
    remainder %= second;
    return remainder;
}


BigInteger& BigInteger::operator++ () {
    *this += 1;
    if (this->GetDigit(0) == 0 && this->GetDigitCount() == 1)
        this->sign = true;
    return *this;
}

BigInteger& BigInteger::operator-- () {
    *this -= 1;
    if (this->GetDigit(0) == 0 && this->GetDigitCount() == 1)
        this->sign = true;
    return *this;
}

const BigInteger BigInteger::operator++ (int) {
    BigInteger copy = *this;
    ++(*this);
    return copy;
}

const BigInteger BigInteger::operator-- (int) {
    BigInteger copy = *this;
    --(*this);
    return copy;
}



Rational::Rational(const int& number) : numerator(number) {}

Rational::Rational(const BigInteger& number) : numerator(number) {}

const BigInteger Rational::GCD(BigInteger first, BigInteger second) const {
    if (first < 0)
        first = -first;
    if (second < 0)
        second = -second;
    while (first != 0 && second != 0) {
        if (first < second)
            second %= first;
        else
            first %= second;
    }
    return (first < second) ? second : first;
}

void Rational::ReductionToIrreducibleForm() {
    BigInteger gcd = GCD(numerator, denominator);
    numerator /= gcd;
    denominator /= gcd;
}

const std::string Rational::toString() const {
    std::string str = numerator.toString();
    if (denominator != 1)
        str += '/' + denominator.toString();
    return str;
}

const std::string Rational::asDecimal(size_t precision = 0) const {
    std::string string;
    BigInteger integer_part = numerator / denominator;
    string += integer_part.toString();
    if (*this < 0 && string[0] != '-')
        string = '-' + string;
    if (precision == 0)
        return string;
    string += '.';
    int ten_pow = 1;
    for (size_t i = 0; i < precision % 9; i++) {
        ten_pow *= 10;
    }
    BigInteger remainder = numerator % denominator;
    if (remainder < 0)
        remainder = -remainder;
    remainder *= ten_pow;
    remainder.DischargeShift(precision / 9);
    std::string str = (remainder / denominator).toString();
    for (size_t i = 0; i < precision - str.size(); i++)
        string += '0';
    string += str;
    return string;
}

Rational::operator double() const {
    return std::stod(numerator.toString())/std::stod(denominator.toString());
}

Rational Rational::operator- () const {
    Rational number = *this;
    number.numerator = -number.numerator;
    return number;
}

bool operator< (const Rational& first, const Rational& second) {
    return (first.numerator * second.denominator < first.denominator * second.numerator);
}

bool operator> (const Rational& first, const Rational& second) {
    return second < first;
}

bool operator<= (const Rational& first, const Rational& second) {
    return !(second < first);
}

bool operator>= (const Rational& first, const Rational& second) {
    return !(first < second);
}

bool operator!= (const Rational& first, const Rational& second) {
    return first < second || second < first;
}

bool operator== (const Rational& first, const Rational& second) {
    return !(first != second);
}

Rational& Rational::operator+= (const Rational& number) {
    this->numerator = this->numerator * number.denominator + this->denominator * number.numerator;
    this->denominator *= number.denominator;
    this->ReductionToIrreducibleForm();
    return *this;
}

Rational& Rational::operator-= (const Rational& number) {
    this->numerator = this->numerator * number.denominator - this->denominator * number.numerator;
    this->denominator *= number.denominator;
    this->ReductionToIrreducibleForm();
    return *this;
}
Rational& Rational::operator*= (const Rational& number) {
    this->numerator *= number.numerator;
    this->denominator *= number.denominator;
    this->ReductionToIrreducibleForm();
    return *this;
}

Rational& Rational::operator/= (const Rational& number) {
    this->numerator *= number.denominator;
    this->denominator *= number.numerator;
    if (this->denominator < 0) {
        this->numerator = -this->numerator;
        this->denominator = -this->denominator;
    }
    this->ReductionToIrreducibleForm();
    return *this;
}

const Rational operator+ (const Rational& first, const Rational& second) {
    Rational sum = first;
    sum += second;
    return sum;
}

const Rational operator- (const Rational& first, const Rational& second) {
    Rational diff = first;
    diff -= second;
    return diff;
}

const Rational operator* (const Rational& first, const Rational& second) {
    Rational mul = first;
    mul *= second;
    return mul;
}

const Rational operator/ (const Rational& first, const Rational& second) {
    Rational quotient = first;
    quotient /= second;
    return quotient;
}
