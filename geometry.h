#include<iostream>
#include<vector>
#include<cmath>
#include<cassert>
#include<cmath>
#include<random>
#include<ctime>

const double EPSILON = 1e-9;
const double INF = 1e60;

bool areEqual(const double first, const double second) {
    return std::abs(first-second) < EPSILON;
}

bool areEqual(const double first, const double second, const double third) {
    return areEqual(first, second) && areEqual(second, third);
}

bool AreProportional(const double first, const double second, const double third, const double fourth) {
    return areEqual(first*fourth, second*third);
}

int signum(double num) {
    if (num > EPSILON)
        return 1;
    if (num < -EPSILON)
        return -1;
    return 0;
}

class Line;

struct Point {
    // точка и вектор вещи схожие, поэтому реализованы в одной сущности
    double x = 0;
    double y = 0;

    Point() = default;
    Point(const double x, const double y) : x(x), y(y) {}

    bool operator== (Point another) const { return areEqual(this->x, another.x) && areEqual(this->y, another.y); }
    bool operator!= (Point another) const { return !(*this == another); }

    double length() const { return sqrt(x*x + y*y); }
    //скалярное произведение
    double operator*(const Point& another) const { return this->x*another.x + this->y*another.y; }
    //умножение на число
    Point operator*(const double num) const { return {x*num, y*num}; }
    Point operator/(const double num) const { return { x/num, y/num}; }
    Point operator+(const Point another) const { return {this->x + another.x, this->y + another.y}; }
    Point operator-(const Point another) const { return {this->x - another.x, this->y - another.y}; }

    void normalize() { *this = *this/length(); }

    void rotate(Point center, double angle);
    void reflex(Point center);
    void reflex(Line axis);
    void scale(Point center, double coefficient);
};

void Point::rotate(const Point center, double angle) {
    angle = angle/180*M_PI;
    double cos = std::cos(angle);
    double sin = std::sin(angle);
    Point diff = *this-center;
    Point new_point(diff.x*cos - diff.y*sin, diff.x*sin + diff.y*cos);
    *this = new_point + center;
}

void Point::reflex(const Point center) {
    *this = center + (center - *this);
}

void Point::scale(const Point center, const double coefficient) {
    *this = center + (*this - center) * coefficient;
}




class Line {
public:
    Line(Point first, Point second);
    Line(double angular_coefficient, double shift);
    Line(Point point, double angular_coefficient);
    Line(double A, double B, double C) : A(A), B(B), C(C) {}

    bool operator== (Line another) const;
    bool operator!= (Line another) const;

    bool IsCorrect() const; // A != 0 or B != 0

    bool hasPoint(Point point) const { return areEqual(0, A*point.x + B*point.y + C); }

    //Ax+By+C = 0
    double A;
    double B;
    double C;
};

bool Line::IsCorrect() const {
    return !areEqual(A, 0, B);
}

Line::Line(const Point first, const Point second) {
    A = first.y-second.y;
    B = second.x-first.x;
    C = first.x*second.y - first.y*second.x;
    assert(IsCorrect());
}

Line::Line(const double angular_coefficient, const double shift) {
    A = angular_coefficient;
    B = -1;
    C = shift;
}

Line::Line(const Point point, const double angular_coefficient) {
    A = angular_coefficient;
    B = -1;
    C = point.y - angular_coefficient*point.x;
}

bool Line::operator== (const Line another) const {
    assert(this->IsCorrect());
    assert(another.IsCorrect());
    return AreProportional(this->A, another.A, this->B, another.B);
}

bool Line::operator!= (const Line another) const {
    return !(*this == another);
}

Line perpendicularLine(Point point, Line line) {
    if (areEqual(line.A, 0))
        return {1, 0, -point.x};
    return {point, line.B/line.A};
}

Point intersectionLines(Line first, Line second) {
    double denominator = first.A*second.B - first.B*second.A;
    if (areEqual(denominator, 0))
        return {INF, INF};
    return {-(first.C*second.B - first.B*second.C) / denominator, -(first.A*second.C - first.C*second.A) / denominator};
}

void Point::reflex(const Line axis) {
    Point projection = intersectionLines(perpendicularLine(*this, axis), axis);
    *this = projection + (projection - *this);
}


double getAngleCosine(Point first, Point second, Point third) {
    Point second_to_first = first-second;
    Point second_to_third = third-second;
    return second_to_first*second_to_third / (second_to_first.length() * second_to_third.length());
}


Line getBisecotor(Point first, Point second, Point third) {
    Point first_vector = first-second;
    Point second_vector = third-second;
    first_vector.normalize();
    second_vector.normalize();
    return {second, second + (first_vector+second_vector)};
}



class Shape {
public:
    virtual ~Shape() = default;

    virtual double perimeter() const = 0;
    virtual double area() const = 0;

    virtual bool operator== (const Shape& another) const = 0;
    bool operator!= (const Shape& another) const { return !(*this == another); }
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;

    virtual void rotate(Point center, double angle) = 0;
    virtual void reflex(Point center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;
};





template<typename T>
void vector_push_back_list(std::vector<T>& vector) {}

template<typename Head, typename... Tail>
void vector_push_back_list(std::vector<Head>& vector, const Head& head, const Tail&... tail) {
    vector.push_back(head);
    vector_push_back_list(vector, tail...);
}

class Polygon : public Shape {
public:
    explicit Polygon(const std::vector<Point>& vertices) : vertices(vertices) {}

    template<typename... List>
    explicit Polygon(const List&... list) { vector_push_back_list(vertices, list...); }

    int verticesCount() const { return vertices.size(); }
    const std::vector<Point>& getVertices() const { return vertices; }

    bool isConvex();

    double perimeter() const override;
    double area() const override;

    bool operator== (const Shape& another) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;

    void rotate(Point center, double angle) override;
    void reflex(Point center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;

protected:
    std::vector<Point> vertices;
    Point getVertex(const int index) const { return vertices.at((index+vertices.size()) % vertices.size()); }
};

int rotateDirection(const Point first, const Point second, const Point third) {
    Point first_to_second = second-first;
    Point second_to_third = third-first;
    return signum(first_to_second.x * second_to_third.y - first_to_second.y * second_to_third.x);
}

bool Polygon::isConvex() {
    int direction = rotateDirection(vertices.back(), vertices[0], vertices[1]); // 1 если по часовой, -1 против часовой
    for (int i = 1; i < verticesCount()-1; i++) {
        if (direction != rotateDirection(vertices[i - 1], vertices[i], vertices[i + 1]))
            return false;
    }
    return direction == rotateDirection(vertices[vertices.size()-2], vertices.back(), vertices[0]);
}

double Polygon::perimeter() const {
    double perimeter = (vertices.back() - vertices[0]).length();
    for (int i = 0; i < verticesCount() - 1; i++)
        perimeter += (vertices[i] - vertices[i+1]).length();
    return perimeter;
}

double Polygon::area() const {
    double area = vertices.back().x * vertices[0].y - vertices[0].x * vertices.back().y;
    for (int i = 0; i < verticesCount()-1; i++) {
        area += vertices[i].x * vertices[i+1].y;
        area -= vertices[i+1].x * vertices[i].y;
    }
    area = std::abs(area)*0.5;
    return area;
}

bool Polygon::operator== (const Shape& another) const {
    auto another_ptr = dynamic_cast<Polygon const*>(&another);
    if (!another_ptr)
        return false;
    if (this->verticesCount() != another_ptr->verticesCount())
        return false;
    const int size = this->verticesCount();

    for (int i = 0; i < size; i++) {
        bool are_equal = true;
        for (int j = 0; j < size; j++) {
            if (this->vertices[j] != another_ptr->vertices[(i + j) % size]) {
                are_equal = false;
                break;
            }
        }
        if (are_equal)
            return true;
    }

    for (int i = 0; i < size; i++) {
        bool are_equal = true;
        for (int j = 0; j < size; j++) {
            if (this->vertices[j] != another_ptr->vertices[(size-1-j+i) % size]) {
                are_equal = false;
                break;
            }
        }
        if (are_equal)
            return true;
    }

    return false;
}

bool Polygon::isSimilarTo(const Shape& another) const {
    auto another_ptr = dynamic_cast<Polygon const*>(&another);
    if (!another_ptr)
        return false;
    if (this->verticesCount() != another_ptr->verticesCount())
        return false;
    const int size = this->verticesCount();

    for (int i = 0; i < size; i++) {
        bool are_similar = true;
        for (int j = 0; j < size; j++) {
            if (!areEqual(getAngleCosine(this->getVertex(j), this->getVertex(j+1), this->getVertex(j+2)),
                          getAngleCosine(another_ptr->getVertex(j+i), another_ptr->getVertex(j+i+1), another_ptr->getVertex(j+i+2)))) {
                are_similar = false;
                break;
            }
        }
        if (are_similar)
            return true;
    }

    for (int i = 0; i < size; i++) {
        bool are_similar = true;
        for (int j = 0; j < size; j++) {
            if (!areEqual(getAngleCosine(this->getVertex(j), this->getVertex(j+1), this->getVertex(j+2)),
                          getAngleCosine(another_ptr->getVertex(-j+i), another_ptr->getVertex(-j+i+1), another_ptr->getVertex(-j+i+2)))) {
                are_similar = false;
                break;
            }
        }
        if (are_similar)
            return true;
    }

    return false;
}

bool Polygon::isCongruentTo(const Shape& another) const {
    return isSimilarTo(another) && areEqual(this->perimeter(), another.perimeter());
}

bool segmentHasPoint(Point first, Point second, Point point) {
    if (first == point || second == point)
        return true;
    if (!Line(first, second).hasPoint(point))
        return false;
    return (point.x > std::min(first.x, second.x) && point.x < std::max(first.x, second.x) &&
            point.y > std::min(first.y, second.y) && point.y < std::max(first.y, second.y));
}

bool beamCrossesSegment(Point beam_start, Point first, Point second) {
    // горизонтальный луч с началом в beam_start направленный в +inf
    // отрезок first <--> second

    assert(!areEqual(beam_start.y, first.y) && !areEqual(beam_start.y, second.y));

    if (!((beam_start.y > first.y) ^ (beam_start.y > second.y)))
        return false;
    double coefficient = std::abs((first.y - beam_start.y)/(second.y - first.y));
    Point point_on_line = first*(1-coefficient) + second*coefficient;
    return point_on_line.x+EPSILON > beam_start.x;
}

bool Polygon::containsPoint(Point point) const {
    for (int i = 0; i < verticesCount(); i++) {
        if (segmentHasPoint(getVertex(i), getVertex(i + 1), point))
            return true;
    }

    // поворот на случайный угол позволяет применить метод подсчёта пересечений рёбер с лучом с началом в данной точке
    Polygon rotated_polygon(*this);
    std::srand(std::time(nullptr));
    rotated_polygon.rotate(point, std::rand()/100.0);
    int intersection_count = 0;
    for  (int i = 0; i < rotated_polygon.verticesCount(); i++) {
        intersection_count += beamCrossesSegment(point, rotated_polygon.getVertex(i), rotated_polygon.getVertex(i+1));
    }
    return intersection_count % 2;
}

void Polygon::rotate(const Point center, double angle) {
    for (Point& vertex : vertices)
        vertex.rotate(center, angle);
}

void Polygon::reflex(const Point center) {
    for (Point& vertex : vertices)
        vertex.reflex(center);
}

void Polygon::reflex(const Line axis) {
    for (Point& vertex : vertices)
        vertex.reflex(axis);
}

void Polygon::scale(Point center, double coefficient) {
    for (Point& vertex : vertices)
        vertex.scale(center, coefficient);
}




class Ellipse : public Shape {
public:
    Ellipse(Point first_focus, Point second_focus, double distance_sum) : first_focus(first_focus), second_focus(second_focus), distance_sum(distance_sum) {}
    std::pair<Point, Point> focuses() const { return {first_focus, second_focus}; }

    double getMinorHalfAxis() const;
    double getMajorHalfAxis() const { return distance_sum*0.5; };
    double getFocalLength() const { return 0.5*(first_focus-second_focus).length(); }

    std::pair<Line, Line> directrices() const;
    double eccentricity() const { return getFocalLength() / getMajorHalfAxis(); }

    double perimeter() const override;
    double area() const override { return M_PI * getMajorHalfAxis() * getMinorHalfAxis(); }

    bool operator== (const Shape& another) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;

    void rotate(Point center, double angle) override;
    void reflex(Point center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;

protected:
    Point first_focus;
    Point second_focus;
    double distance_sum;

};

double Ellipse::getMinorHalfAxis() const {
    double a = getMajorHalfAxis();
    double c = getFocalLength();
    return sqrt(a*a-c*c);
}

std::pair<Line, Line> Ellipse::directrices() const {
    double b = getMinorHalfAxis();
    double distance = b*b / (eccentricity() * getMajorHalfAxis()); // расстояние между фокусом и директрисой

    Point vector = first_focus-second_focus; // second to first
    vector.normalize();

    Line majorAxis(first_focus, second_focus);

    return {perpendicularLine(first_focus + vector*distance, majorAxis),
            perpendicularLine(second_focus - vector*distance, majorAxis)};
}


double Ellipse::perimeter() const {
    double a = getMajorHalfAxis();
    double b = getMinorHalfAxis();
    return M_PI*(3*(a+b) - sqrt((3*a+b) * (a+3*b)));
}

bool Ellipse::operator==(const Shape& another) const {
    auto another_ptr = dynamic_cast<Ellipse const*>(&another);
    if (!another_ptr)
        return false;
    if (!areEqual(this->distance_sum, another_ptr->distance_sum))
        return false;
    return (this->first_focus == another_ptr->first_focus && this->second_focus == another_ptr->second_focus) ||
           (this->first_focus == another_ptr->second_focus && this->second_focus == another_ptr->first_focus);
}

bool Ellipse::isCongruentTo(const Shape& another) const {
    auto another_ptr = dynamic_cast<Ellipse const*>(&another);
    if (!another_ptr)
        return false;
    if (!areEqual(this->distance_sum, another_ptr->distance_sum))
        return false;
    return areEqual((this->first_focus - this->second_focus).length(),
                    (another_ptr->first_focus - another_ptr->second_focus).length());
}

bool Ellipse::isSimilarTo(const Shape &another) const {
    auto another_ptr = dynamic_cast<Ellipse const*>(&another);
    if (!another_ptr)
        return false;
    return areEqual(this->distance_sum / another_ptr->distance_sum, (this->first_focus - this->second_focus).length() / (another_ptr->first_focus - another_ptr->second_focus).length());
}

bool Ellipse::containsPoint(Point point) const {
    return (point - this->first_focus).length() + (point - this->second_focus).length() < distance_sum + EPSILON;
}

void Ellipse::rotate(Point center, double angle) {
    first_focus.rotate(center, angle);
    second_focus.rotate(center, angle);
}

void Ellipse::reflex(Point center) {
    first_focus.reflex(center);
    second_focus.reflex(center);
}

void Ellipse::reflex(Line axis) {
    first_focus.reflex(axis);
    second_focus.reflex(axis);
}

void Ellipse::scale(Point center, double coefficient) {
    first_focus.scale(center, coefficient);
    second_focus.scale(center, coefficient);
    distance_sum *= std::abs(coefficient);
}




class Circle : public Ellipse {
public:
    Circle(Point center, double radius) : Ellipse(center, center, 2*radius) {}

    bool isSimilarTo(const Shape& another) const override;

    Point center() { return first_focus; }
    double radius() const { return getMajorHalfAxis(); }
    double perimeter() const override { return 2*M_PI*radius(); }
};

bool Circle::isSimilarTo(const Shape& another) const {
    auto another_ptr = dynamic_cast<Circle const*>(&another);
    return another_ptr != nullptr;
}


class Rectangle : public Polygon {
public:
    template<typename... List>
    explicit Rectangle(const List&... list) : Polygon(list...) {}

    explicit Rectangle(const std::vector<Point>& vertices) : Polygon(vertices) {}

    Rectangle(Point first, Point second, double coefficient);
    Rectangle(Point first, Point second, int coefficient) : Rectangle(first, second, static_cast<double>(coefficient)) {}

    Point center() const { return (vertices[0]+vertices[2])*0.5; }
    std::pair<Line, Line> diagonals() const { return {{vertices[0], vertices[2]}, {vertices[1], vertices[3]}}; }
};

Rectangle::Rectangle(const Point first, const Point second, const double coefficient) {
    vertices.resize(4);
    vertices[0] = first;
    vertices[2] = second;
    Point vector = second-first;
    double angle = atan(std::max(coefficient,1/coefficient));
    vector.rotate({0, 0}, angle/M_PI*180);
    double smaller_edge_length = vector.length()*std::min(cos(angle), sin(angle));
    vector = vector*(smaller_edge_length/(vector.length()));
    vertices[1] = vertices[0] + vector;
    vertices[3] = vertices[2] - vector;
}




class Square : public Rectangle {
public:
    // не знаю как лучше делать: так, или как в Rectangle
    Square(Point first, Point second) : Rectangle(first, second, 1.0) {}
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

Circle Square::circumscribedCircle() const {
    Point circle_center = center();
    return {circle_center, (vertices[0]-circle_center).length()};
}

Circle Square::inscribedCircle() const {
    Point circle_center = center();
    Point edge_mid = (vertices[0] + vertices[1]) * 0.5;
    return {circle_center, (edge_mid-circle_center).length()};
}




class Triangle : public Polygon {
public:
    Triangle(Point first, Point second, Point third) : Polygon(first, second, third) {}

    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
    Point centroid() const;
    Point orthocenter() const;
    Line EulerLine() const;
    Circle ninePointsCircle() const;
};

Circle Triangle::circumscribedCircle() const {
    Point center = intersectionLines(perpendicularLine((vertices[0]+vertices[1]) * 0.5, {vertices[0], vertices[1]}),
                                     perpendicularLine((vertices[1]+vertices[2]) * 0.5, {vertices[1], vertices[2]}));
    return {center, (vertices[0] - center).length()};
}

Circle Triangle::inscribedCircle() const {
    Point center = intersectionLines(getBisecotor(vertices[0], vertices[1], vertices[2]),
                                     getBisecotor(vertices[1], vertices[2], vertices[0]));
    return {center, 2*area()/perimeter()};
}

Point Triangle::centroid() const {
    Point center = intersectionLines({vertices[0], (vertices[1] + vertices[2]) * 0.5},
                                     {vertices[1], (vertices[0] + vertices[2]) * 0.5});
    assert(center != Point(INF, INF));
    return center;
}

Point Triangle::orthocenter() const {
    Point center = intersectionLines(perpendicularLine(vertices[0], {vertices[1], vertices[2]}),
                                     perpendicularLine(vertices[1], {vertices[0], vertices[2]}));
    assert(center != Point(INF, INF));
    return center;
}

Line Triangle::EulerLine() const {
    return {orthocenter(), centroid()};
}

Circle Triangle::ninePointsCircle() const {
    return Triangle((vertices[0]+vertices[1])*0.5, (vertices[1]+vertices[2])*0.5, (vertices[2]+vertices[0])*0.5).circumscribedCircle();
}
