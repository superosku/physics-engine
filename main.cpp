#include <iostream>
#include <vector>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <stdlib.h>
#include <math.h>
#include <cmath>


float zoom = 5.0;


class Vec2f {
protected:
public:
    float x;
    float y;

    Vec2f(float x, float y) {
        this->x = x;
        this->y = y;
    }

    float distance(Vec2f *other) {
        float diffx = this->x - other->x;
        float diffy = this->y - other->y;
        return sqrt(diffx * diffx + diffy * diffy);
    }

    float length() {
        return sqrt(this->x * this->x + this->y * this->y);
    }

    Vec2f minus(Vec2f *other) {
        return Vec2f(this->x - other->x, this->y - other->y);
    }

    Vec2f added(Vec2f *other) {
        return Vec2f(this->x + other->x, this->y + other->y);
    }

    Vec2f multiply(float value) {
        return Vec2f(this->x * value, this->y * value);
    }

    Vec2f normalized() {
        float dist = sqrt(this->x * this->x + this->y * this->y);
        return Vec2f(this->x / dist, this->y / dist);
    }
};

class Circle {
protected:
public:
    Vec2f center;
    float radius;

    Circle(Vec2f center, float radius) : center(0.0, 0.0) {
        this->radius = radius;
        this->center = center;
    }
};

class BoundingBox {
public:
    float xMin;
    float xMax;
    float yMin;
    float yMax;

    bool isInside(Vec2f point) {
        return (
            point.x >= this->xMin && point.x <= this->xMax &&
            point.y >= this->yMin && point.y <= this->yMax
        );
    }
};

bool intersectsVertical(
    Vec2f point1,
    Vec2f point2,
    float x1,
    float x2,
    float height
) {
    // Two vertical lines do not intersect
    if (point1.y == point2.y) {
        return false;
    }
    // Height is out of bounds
    if (
        height >= std::max(point1.y, point2.y) ||
        height <= std::min(point1.y, point2.y)
    ) {
        return false;
    }
    // X direction is out of bounds
    if (
        std::min(x1, x2) >= std::max(point1.x, point2.x) ||
        std::max(x1, x2) <= std::min(point1.x, point2.x)
    ) {
        return false;
    }

    // Inside bounds and does intersect
    float miny = std::min(point1.y, point2.y);
    float minx = point1.y < point2.y ? point1.x : point2.x;
    float minratio = std::abs(miny - height);
    float maxy = std::max(point1.y, point2.y);
    float maxx = point1.y > point2.y ? point1.x : point2.x;
    float maxratio = std::abs(maxy - height);
    float xPos = minx + (maxx - minx) * (minratio / (maxratio + minratio));

    return std::min(x1, x2) < xPos && std::max(x1, x2) > xPos;
}

struct PointInfo {
    Vec2f point;
    bool isOnePoint;
    float d1;
    float d2;
};

PointInfo closestPointToLine(
    Vec2f *p1,
    Vec2f *p2,
    Vec2f *point
) {
    if (p1->x == p2->x) {
        // TODO
        return PointInfo{Vec2f(p1->x, point->y), false, std::abs(point->y - p1->y), std::abs(point->y - p2->y)};
    }
    if (p1->y == p2->y) {
        // TODO
        return PointInfo{Vec2f(point->x, p1->y), false, std::abs(point->x - p1->x), std::abs(point->x - p2->x)};
    }
    Vec2f dirPoints = p1->minus(p2);
    Vec2f dirMove = Vec2f(dirPoints.y, -dirPoints.x).normalized();

    // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
    float distance = (
        (
            (p2->x - p1->x) * (p1->y - point->y) -
            (p1->x - point->x) * (p2->y - p1->y)
        ) /
        sqrt(
            (p2->x - p1->x) * (p2->x - p1->x) +
            (p2->y - p1->y) * (p2->y - p1->y)
        )
    );

    if ( sqrt(
            (p2->x - p1->x) * (p2->x - p1->x) +
            (p2->y - p1->y) * (p2->y - p1->y)
        ) == 0.0
    ) {
        std::cout << "HMMMMM" << std::endl;
    }

    Vec2f move = dirMove.multiply(distance);

    Vec2f possiblePoint = point->added(&move);

    // if point is not in line
    if (
        possiblePoint.x < std::min(p1->x, p2->x) ||
        possiblePoint.x > std::max(p1->x, p2->x) ||
        possiblePoint.y < std::min(p1->y, p2->y) ||
        possiblePoint.y > std::max(p1->y, p2->y)
    ) {
        // No need to check if we should return this from p1 or p2 since
        // The the point is in two lines and always one of them will be the
        // right one.. :)
        return PointInfo{Vec2f(p1->x, p1->y), true, 1.0, 0.0};
    }

    // This might prevent bodies from sliding into each other
//    if ( possiblePoint.distance(p1) < 0.5 ) {
//        return PointInfo{ *p1, false, 1.0, 0.0 };
//    }
//    if ( possiblePoint.distance(p2) < 0.5 ) {
//        return PointInfo{ *p2, false, 0.0, 1.0 };
//    }

    return PointInfo{
        possiblePoint,
        false,
        possiblePoint.distance(p1),
        possiblePoint.distance(p2)
    };
}

class VerletPoint {
public:
    Vec2f pos;
    Vec2f prevPos;
    Vec2f acceleration;
    float weight;

    VerletPoint(Vec2f position)
    : pos{position},
      prevPos{position},
      acceleration{0.0, 0.0},
      weight{1}
    { }

    void update(float subStep) {
        Vec2f displacement = this->pos.minus(&this->prevPos);
        this->prevPos = this->pos;
//        Vec2f accSS = this->acceleration.multiply(subStep * subStep);
        Vec2f accSS = this->acceleration.multiply(subStep);
        if (displacement.length() > 1.1) {
            std::cout << "WTWTWTWTF1 " << displacement.length() << std::endl;
        }
        if (accSS.length() > 1.1) {
            std::cout << "WTWTWTWTF2 " << accSS.length() << std::endl;
        }
        this->pos = this->pos
            .added(&displacement)
            .added(&accSS);
        this->acceleration = Vec2f(0.0, 0.0);
    }

    void accelerate(Vec2f acceleration) {
        this->acceleration = this->acceleration.added(&acceleration);
    }

    void addVelocity(Vec2f velocity) {
        this->prevPos = this->prevPos.minus(&velocity);
    }

    void setVelocity(Vec2f velocity) {
        this->prevPos = velocity;
    }
};


class Spring {
public:
    int i1;
    int i2;
    float length;

    Spring(int i1, int i2, float length) {
        this->i1 = i1;
        this->i2 = i2;
        this->length = length;
    }
};


struct ClosestPointInfo {
    Vec2f point;
    int i;
    int j;
    float distance;
    float d1;
    float d2;
};


struct SpringInfo {
    int i;
    int j;
};


struct IntersectionStatus {
    bool intersects;
    float x;
    float y;
};


IntersectionStatus linesIntersect(Vec2f& point11, Vec2f& point12, Vec2f& point21, Vec2f& point22) {
    float x1 = point11.x;
    float y1 = point11.y;
    float x2 = point12.x;
    float y2 = point12.y;
    float x3 = point21.x;
    float y3 = point21.y;
    float x4 = point22.x;
    float y4 = point22.y;

    float den = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
    if (den == 0.0) {
        // The two lines are parallel and do not intersect
        return IntersectionStatus{false, 0.0, 0.0};
//        return false;
    }

    float ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / den;
    float ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / den;

    if (ua < 0.0 || ua > 1.0 || ub < 0.0 || ub > 1.0) {
        // The intersection point is not within the line segments
        return IntersectionStatus{false, 0.0, 0.0};
//        return false;
    }

    // Calculate the intersection point
    float x = x1 + ua * (x2 - x1);
    float y = y1 + ua * (y2 - y1);

    // The lines intersect at (x, y)
    return IntersectionStatus{true, x, y};
//    return true;
}

float angleBetweenLines(Vec2f& line1point1, Vec2f& line1point2, Vec2f& line2point1, Vec2f& line2point2) {
    Vec2f line1 = {line1point2.x - line1point1.x, line1point2.y - line1point1.y};
    Vec2f line2 = {line2point2.x - line2point1.x, line2point2.y - line2point1.y};

    float dotProduct = line1.x * line2.x + line1.y * line2.y;
    float magLine1 = std::sqrt(line1.x * line1.x + line1.y * line1.y);
    float magLine2 = std::sqrt(line2.x * line2.x + line2.y * line2.y);
    float cosTheta = dotProduct / (magLine1 * magLine2);
    float angle = std::acos(cosTheta);

    return angle;
}

class Body {
protected:
public:
    std::vector<VerletPoint> points;
    std::vector<Spring> springs;
    bool isStatic;

    Body(std::vector<Vec2f> points, bool isStatic) {
        // Is static
        this->isStatic = isStatic;
        // Points
        for (auto &point: points) {
            this->points.push_back(VerletPoint(point));
        }
        // Springs
        if (isStatic) {
            return;
        }
        for (int i = 0; i < points.size(); i ++) {
            Vec2f point1 = points[i];
            int j = (i + 1) % points.size();
            Vec2f point2 = points[j];
            this->springs.push_back(Spring(i, j, point1.distance(&point2)));

            int k = (i + 2) % points.size();
            Vec2f point3 = points[k];
            this->springs.push_back(Spring(i, k, point1.distance(&point3)));
        }
    }

    Body(std::vector<Vec2f> points, std::vector<SpringInfo> springInfos) {
        // Is static
        this->isStatic = false;
        // Points
        for (auto &point: points) {
            this->points.push_back(VerletPoint(point));
        }
        // Springs
        for (auto &springInfo: springInfos) {
            Vec2f point1 = points[springInfo.i];
            Vec2f point2 = points[springInfo.j];
            this->springs.push_back(Spring(springInfo.i, springInfo.j, point1.distance(&point2)));
        }
    }

    BoundingBox getBoundingBox() {
        BoundingBox box;
        box.xMax = this->points[0].pos.x;
        box.xMin = this->points[0].pos.x;
        box.yMax = this->points[0].pos.y;
        box.yMin = this->points[0].pos.y;
        for (auto &point : this->points) {
            box.xMax = std::max(box.xMax, point.pos.x);
            box.xMin = std::min(box.xMin, point.pos.x);
            box.yMax = std::max(box.yMax, point.pos.y);
            box.yMin = std::min(box.yMin, point.pos.y);
        }
        return box;
    }

    bool isInside(Vec2f point) {
        int intersectionCount = 0;

        BoundingBox boundingBox = this->getBoundingBox();

        // Quick hack that might fix issue with explosions
        if (point.x < boundingBox.xMin || point.x > boundingBox.xMax) {
            return false;
        }
        if (point.y < boundingBox.yMin || point.y > boundingBox.yMax) {
            return false;
        }

        for (int i = 0; i < this->points.size(); i++) {
            Vec2f point1 = this->points[i].pos;
            Vec2f point2 = this->points[(i + 1) % this->points.size()].pos;
            if (intersectsVertical(
                point1,
                point2,
                point.x,
                boundingBox.xMax + 1.0,
                point.y
            )) {
                intersectionCount += 1;
            }
        }

        return intersectionCount % 2 == 1;
    }

    ClosestPointInfo closestPoint(Vec2f &point) {
        struct ClosestPointInfo best{
            Vec2f(0.0, 0.0),
            0,
            0,
            1000000000000.0,
            0.5,
            0.5
        };

        for (int i = 0; i < this->points.size(); i ++) {
            Vec2f point1 = this->points[i].pos;
            int j = (i + 1) % this->points.size();
            Vec2f point2 = this->points[j].pos;

            PointInfo pointCandidate = closestPointToLine(&point1, &point2, &point);
            float candidateDistance = pointCandidate.point.distance(&point);

            if (candidateDistance < best.distance) {
                best.distance = candidateDistance;
                best.point = pointCandidate.point;
                best.i = i;
                best.j = j;
                best.d1 = pointCandidate.d1;
                best.d2 = pointCandidate.d2;
            }
        }

        return best;
    }
};

class String {
    std::vector<Vec2f> Vec2fs;
public:
    String(std::vector<Vec2f> Vec2fs) {
        this->Vec2fs = Vec2fs;
    }
};

class Simulation {

private:
public:
//    std::vector<Circle> circles;
    std::vector<Body> bodies;

    Simulation() {
    }

//    void addCircle(Circle circle) {
//        this->circles.push_back(circle);
//    }

    void addBody(Body body) {
        this->bodies.push_back(body);
    }

    void applyGravity(float subStep) {
        for (int i = 0; i < this->bodies.size(); i ++) {
            Body& body1 = this->bodies[i];
            if (body1.isStatic) {
                continue;
            }
            for (auto &point: body1.points) {
                point.accelerate(Vec2f(0.0, 0.05 * subStep));
            }
        }
    }

    void applySprings(float subStep) {
        for (int i = 0; i < this->bodies.size(); i ++) {
            Body& body1 = this->bodies[i];
            if (body1.isStatic) {
                continue;
            }

            for (auto &spring: body1.springs) {
                auto *point1 = &body1.points[spring.i1];
                auto *point2 = &body1.points[spring.i2];
                float distance = point1->pos.distance(&point2->pos);
                float forceScalar = (distance - spring.length);
                Vec2f forceDirection = point1->pos.minus(&point2->pos).normalized();

                // Accelerate
//                point1->accelerate(forceDirection.multiply(forceScalar * -0.5 * subStep * 0.2));
//                point2->accelerate(forceDirection.multiply(forceScalar * 0.5 * subStep * 0.2));
                // Move
                float moveFactor = 0.06 * subStep; // For springiness
                Vec2f v1 = forceDirection.multiply(forceScalar * -0.5 * 1.0 * 0.75 * moveFactor);
                point1->pos = point1->pos.added(&v1);
                Vec2f v2 = forceDirection.multiply(forceScalar * 0.5 * 1.0 * 0.75 * moveFactor);
                point2->pos = point2->pos.added(&v2);

                if (v1.length() > 1.00 || v2.length() > 1.00) {
                    std::cout << "WTF " << i << " " << v1.length() << " " << v2.length() << std::endl;
                }
            }
        }
    }
    void applyCollisions(float subStep) {
        for (int i = 0; i < this->bodies.size(); i ++) {
            Body& body1 = this->bodies[i];
            if (body1.isStatic) {
                continue;
            }

            for (int j = 0; j < this->bodies.size(); j ++) {
                if (i == j) {
                    continue;
                }
                Body& body2 = this->bodies[j];

                for (auto &point: body1.points) {
                    bool isInside = body2.isInside(point.pos);
                    if (isInside) {
//                        if (i > 2 && j >> 2) {
//                            std::cout << "is inside " << i << " " << j << std::endl;
//
//                            std::cout << "Object1: ";
//                            for (auto &ppp: body1.points) {
//                                std::cout << "(" << ppp.pos.x << "," << ppp.pos.y << ") ";
//                            }
//                            std::cout << std::endl;
//
//                            std::cout << "Object1 point: ";
//                            std::cout << "(" << point.pos.x << "," << point.pos.y << ") ";
//                            std::cout << std::endl;
//
//                            std::cout << "Object2: ";
//                            for (auto &ppp: body2.points) {
//                                std::cout << "(" << ppp.pos.x << "," << ppp.pos.y << ") ";
//                            }
//                            std::cout << std::endl;
//                        }
                        ClosestPointInfo pointInfo = body2.closestPoint(point.pos);

                        Vec2f direction = pointInfo.point.minus(&point.pos).multiply(0.5 * 0.75);

                        if (direction.length() > 1.0) {
                            std::cout << "direction len " << direction.length() << std::endl;
                        }

                        // Move the inside point
                        point.pos = point.pos.added(&direction);
//                        point.addVelocity(direction);
//                        point.setVelocity(pointInfo.point);

                        // Move the two points in the line that was overlapped
                        if (body2.isStatic) {
                            continue;
                        }
                        auto &point1 = body2.points[pointInfo.i];
                        auto &point2 = body2.points[pointInfo.j];
                        if ((pointInfo.d1 + pointInfo.d2) == 0) {
                            std::cout << "HMMMMMMMM" << std::endl;
                            continue;
                        }
                        Vec2f od1 = direction.multiply(-1.0 * (pointInfo.d2 / (pointInfo.d1 + pointInfo.d2)));
                        Vec2f od2 = direction.multiply(-1.0 * (pointInfo.d1 / (pointInfo.d1 + pointInfo.d2)));
                        point1.pos = point1.pos.added(&od1);
                        point2.pos = point2.pos.added(&od2);
//                        Vec2f otherDirection = direction.multiply(-1.0);
//                        point1.addVelocity(otherDirection);
//                        point2.addVelocity(otherDirection);
                    }
                }
            }

            for (int j = 0; j < this->bodies.size(); j ++) {
                if (j == i) {
                    continue;
                }
                Body& body2 = this->bodies[j];

                if (body1.isStatic && body2.isStatic) {
                    continue;
                }

                for (int body1point1index = 0; body1point1index < body1.points.size(); body1point1index ++) {
                    int body1point2index = (body1point1index + 1) % body1.points.size();
                    for (int body2point1index = 0; body2point1index < body2.points.size(); body2point1index ++) {
                        int body2point2index = (body2point1index + 1) % body2.points.size();
                        VerletPoint& point11 = body1.points[body1point1index];
                        VerletPoint& point12 = body1.points[body1point2index];
                        VerletPoint& point21 = body2.points[body2point1index];
                        VerletPoint& point22 = body2.points[body2point2index];

                        IntersectionStatus intersectStatus = linesIntersect(
                            point11.pos, point12.pos, point21.pos, point22.pos
                        );

                        float angle = angleBetweenLines(point11.pos, point12.pos, point21.pos, point22.pos);
//                        std::cout << "angle: " << angle << std::endl;
                        float angleCutoff = 5.0;

                        if (intersectStatus.intersects && angle > angleCutoff && angle < 6.28318531 - angleCutoff) {
                            std::cout << "There is an intersection " << i << " " << j << " " << angle << std::endl;

                            Vec2f intersectionPoint = Vec2f(intersectStatus.x, intersectStatus.y);

                            float distToP11 = intersectionPoint.distance(&point11.pos);
                            float distToP12 = intersectionPoint.distance(&point12.pos);
                            float distToP21 = intersectionPoint.distance(&point21.pos);
                            float distToP22 = intersectionPoint.distance(&point22.pos);

                            if (!body1.isStatic && distToP11 < 1.1) {
//                            if (distToP11 < distToP12 && distToP11 < distToP21 && distToP11 < distToP22 && !body1.isStatic && distToP11 < 0.5) {
                                point11.pos.x = intersectionPoint.x;
                                point11.pos.y = intersectionPoint.y;
                            }
//                            if (distToP12 < distToP11 && distToP12 < distToP21 && distToP12 < distToP22 && !body1.isStatic && distToP12 < 0.5) {
                            if (!body1.isStatic && distToP12 < 1.1) {
                                point12.pos.x = intersectionPoint.x;
                                point12.pos.y = intersectionPoint.y;
                            }
//                            if (distToP21 < distToP12 && distToP21 < distToP11 && distToP21 < distToP22 && !body2.isStatic && distToP21 < 0.5) {
                            if (!body2.isStatic && distToP21 < 1.1) {
                                point21.pos.x = intersectionPoint.x;
                                point21.pos.y = intersectionPoint.y;
                            }
//                            if (distToP22 < distToP12 && distToP22 < distToP21 && distToP22 < distToP11 && !body2.isStatic && distToP22 < 0.5) {
                            if (!body2.isStatic && distToP22 < 1.1) {
                                point22.pos.x = intersectionPoint.x;
                                point22.pos.y = intersectionPoint.y;
                            }
                        }
                    }
                }
            }
        }
    }
    void updateObjects(float subStep) {
        for (int i = 0; i < this->bodies.size(); i ++) {
            Body& body1 = this->bodies[i];
            if (body1.isStatic) {
                continue;
            }

            for (auto &point: body1.points) {
                point.update(subStep);
            }
        }
    }

    void simulate() {
        int subSteps = 4;
        float step = 1.0 / subSteps;
        for (int i = 0; i < subSteps ; i ++) {
            this->applyGravity(step);
            this->applySprings(step);
            this->applyCollisions(step);
            this->updateObjects(step);
        }

//        for (auto & body: this->bodies) {
//            auto boundingBox = body.getBoundingBox();
//            std::cout << boundingBox.xMin << " " << boundingBox.xMax << std::endl;
//        }
    }

    void draw(sf::RenderWindow *window) {

//        for (auto &circle: this->circles) {
//            sf::CircleShape shape(circle.radius * zoom);
//            shape.setPosition(circle.center.x * zoom, circle.center.y * zoom);
//            shape.setFillColor(sf::Color::Green);
//
//            window->draw(shape);
//        }

        for (auto &body: this->bodies) {
            // Outline
            sf::ConvexShape convex;
            convex.setPointCount(body.points.size());
            convex.setOutlineThickness(4);
            convex.setOutlineColor(sf::Color(128, 128, 128));
            convex.setFillColor(sf::Color::Transparent);
            for (int i = 0; i < body.points.size(); i ++) {
                Vec2f point = body.points[i].pos;

                convex.setPoint(
                    i,
                    sf::Vector2f(
                        body.points[i].pos.x * zoom,
                        body.points[i].pos.y * zoom
                    )
                );
            }
            window->draw(convex);

            // Corner points
            for (int i = 0; i < body.points.size(); i ++) {
                Vec2f point = body.points[i].pos;

                sf::CircleShape shape(1.5 * zoom);
                shape.setPosition(
                    point.x * zoom - 1.5 * zoom,
                    point.y * zoom - 1.5 * zoom
                );
                shape.setFillColor(sf::Color::Green);

                window->draw(shape);
            }
        }
    }
};


void addRectangle(Simulation &simulation, float x, float y, float size) {
    std::vector<Vec2f> bodyPoints {
        Vec2f(x + size, y + size),
        Vec2f(x + 0.0, y + size),
        Vec2f(x + 0.0, y + 0.0),
        Vec2f(x + size, y + 0.0),
    };
    simulation.addBody(Body(bodyPoints, false));
}

void addTriangle(Simulation &simulation, float x, float y, float size) {
    std::vector<Vec2f> bodyPoints {
        Vec2f(x + size, y + size),
        Vec2f(x + 0.0, y + size),
        Vec2f(x + size / 2.0, y + 0.0),
    };
    simulation.addBody(Body(bodyPoints, false));
}


void addRectangle2(Simulation &simulation, float x, float y, float size) {
    std::vector<Vec2f> bodyPoints {
        Vec2f(x + size * 1.0, y + size * 1.0),
        Vec2f(x + size * 0.5, y + size * 1.0),
        Vec2f(x + size * 0.0, y + size * 1.0),
        Vec2f(x + size * 0.0, y + size * 0.5),
        Vec2f(x + size * 0.0, y + size * 0.0),
        Vec2f(x + size * 0.5, y + size * 0.0),
        Vec2f(x + size * 1.0, y + size * 0.0),
        Vec2f(x + size * 1.0, y + size * 0.5),
    };
    std::vector<SpringInfo> springInfos;
    for (int i = 0; i < 8; i ++) {
        int j = (i + 1) % 8;
//        int k = (i + 2) % 8;
        int l = (i + 3) % 8;
        int m = (i + 4) % 8;
//        int n = (i + 5) % 8;
        springInfos.push_back(SpringInfo{i, j});
//        springInfos.push_back(SpringInfo{i, k});
        springInfos.push_back(SpringInfo{i, l});
        springInfos.push_back(SpringInfo{i, m});
//        springInfos.push_back(SpringInfo{i, n});
    }
    simulation.addBody(Body(bodyPoints, springInfos));
}


int main(int argc, char *argv[]) {

    sf::RenderWindow window(sf::VideoMode(2000, 1600), "SFML works!");

    Simulation simulation;

    {
        std::vector<Vec2f> bodyPoints {
            Vec2f(10.0, 250.0),
            Vec2f(10.0, 310.0),
            Vec2f(390.0, 310.0),
            Vec2f(390.0, 250.0),
        };
        simulation.addBody(Body(bodyPoints, true));
    }
    {
        std::vector<Vec2f> bodyPoints {
            Vec2f(10.0, 250.0),
            Vec2f(10.0, 110.0),
            Vec2f(30.0, 110.0),
            Vec2f(30.0, 250.0),
        };
        simulation.addBody(Body(bodyPoints, true));
    }
    {
        std::vector<Vec2f> bodyPoints {
            Vec2f(370.0, 250.0),
            Vec2f(370.0, 110.0),
            Vec2f(390.0, 110.0),
            Vec2f(390.0, 250.0),
        };
        simulation.addBody(Body(bodyPoints, true));
    }

//    addTriangle(simulation, 50.0, 195.0, 50.0);
//    addTriangle(simulation, 150.0, 195.0, 50.0);
//    addRectangle(simulation, 190.0, 210.0, 10.0);
//    addRectangle(simulation, 200.0, 190.0, 10.0);

    for (int i = 0; i < 6; i ++) {
        for (int j = 0; j < 6; j ++) {
//            if ((i + j) % 2 == 0) {
                addRectangle(
                    simulation,
                    60.0 + i * 40.0 + (j % 2 == 1 ? 15.0 : 0.0),
                    50.0 + j * 50.0 - 100.0,
                    20.0 + i * 2.0 + j * 3.0
                );
//            } else {
//                addTriangle(
//                    simulation,
//                    60.0 + i * 40.0 + (j % 2 == 1 ? 15.0 : 0.0),
//                    50.0 + j * 40.0,
//                    20.0 + i * 2.0 + j * 3.0
//                );
//            }
        }
    }

//    addRectangle(
//        simulation,
//        60.0,
//        190.0,
//        50.0
//    );
//    addRectangle(
//        simulation,
//        65.0,
//        100.0,
//        50.0
//    );
//    addRectangle(
//        simulation,
//        265.0,
//        100.0,
//        50.0
//    );

//    simulation.addCircle(Circle(Vec2f(1.0, 1.0), 1.0));
//    simulation.addCircle(Circle(Vec2f(5.0, 5.0), 2.0));
//    {
//        std::vector<Vec2f> bodyPoints {
//            Vec2f(120.0, 20.0),
//            Vec2f(180.0, 40.0),
////            Vec2f(200.0, 20.0),
////            Vec2f(160.0, 50.0),
////            Vec2f(180.0, 80.0),
//            Vec2f(140.0, 90.0),
////            Vec2f(120.0, 80.0),
//        };
//        simulation.addBody(Body(bodyPoints, false));
//    }
//    {
//        std::vector<Vec2f> bodyPoints {
//            Vec2f(300.0, 40.0),
//            Vec2f(300.0, 0.0),
//            Vec2f(340.0, 40.0),
//        };
//        simulation.addBody(Body(bodyPoints, false));
//    }

    sf::Clock clock;
    sf::Time timeSinceLastUpdate = sf::Time::Zero;

    sf::Time TimePerFrame = sf::seconds(1.f / 60.f);
//    sf::Time TimePerFrame = sf::seconds(1.f / 2.f);

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        timeSinceLastUpdate += clock.restart();
        while (timeSinceLastUpdate > TimePerFrame)
        {
            timeSinceLastUpdate -= TimePerFrame;

            window.clear();

            simulation.simulate();

            simulation.draw(&window);

            // Draw mouse
            sf::Vector2i mousePosition = sf::Mouse::getPosition(window);
            Vec2f simulationMousePos = Vec2f(mousePosition.x / zoom, mousePosition.y / zoom);
            bool isInside = simulation.bodies[0].isInside(simulationMousePos);
            float radius = isInside ? 20.0 : 10.0;
            sf::CircleShape shape(radius);
            shape.setPosition(
                mousePosition.x - radius,
                mousePosition.y - radius
            );
            shape.setFillColor(isInside ? sf::Color::Green: sf::Color::Red);
            window.draw(shape);

            if (isInside) {
                ClosestPointInfo closeInfo = simulation.bodies[0].closestPoint(simulationMousePos);

                sf::CircleShape shape(10.0);
                shape.setPosition(
                    closeInfo.point.x * zoom - 10.0,
                    closeInfo.point.y * zoom - 10.0
                );
                shape.setFillColor(isInside ? sf::Color::Green: sf::Color::Yellow);
                window.draw(shape);
            }

            for (int i = 0; i < simulation.bodies.size(); i ++) {
                Body& body1 = simulation.bodies[i];

                for (int body1point1index = 0; body1point1index < body1.points.size(); body1point1index ++) {
                    int body1point2index = (body1point1index + 1) % body1.points.size();

                    VerletPoint& point11 = body1.points[body1point1index];
                    VerletPoint& point12 = body1.points[body1point2index];

                    auto origo = Vec2f(0.0, 0.0);
                    IntersectionStatus intersectStatus = linesIntersect(
                        point11.pos, point12.pos, simulationMousePos, origo
                    );

                    if (intersectStatus.intersects) {
                        sf::CircleShape shape(10.0);
                        shape.setPosition(
                            intersectStatus.x * zoom - 10.0,
                            intersectStatus.y * zoom - 10.0
                        );
                        shape.setFillColor(sf::Color(255, 0, 255));
                        window.draw(shape);
                    }
                }
            }

    //        // Draw
    //        for (int i = 0; i < simulation.bodies[0].points.size(); i ++) {
    //            Vec2f point1 = simulation.bodies[0].points[i];
    //            Vec2f point2 = simulation.bodies[0].points[(i + 1) % simulation.bodies[0].points.size()];
    //
    //            Vec2f closestPoint = closestPointToLine(&point1, &point2, &simulationMousePos);
    //
    //            sf::CircleShape shape(10.0);
    //            shape.setPosition(
    //                closestPoint.x * zoom - 10.0,
    //                closestPoint.y * zoom - 10.0
    //            );
    //            shape.setFillColor(isInside ? sf::Color::Green: sf::Color::Yellow);
    //            window.draw(shape);
    //        }

            window.display();
        }
    }

    return 0;

}
