// レンダリング

#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>
#include <map>
#include <chrono>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>

using namespace Eigen;

struct Halfedge {
    int idx;
    int h_opp;

    int face(){
        return int(idx/3);
    }
    int h_next(){
        if(idx % 3 == 2) return idx-2;
        else return idx+1;
    }
    int h_prev(){
        if(idx % 3 == 0) return idx+2;
        else return idx-1;
    }
    int v_src(std::vector<Vector3i> &faces){
        return faces[face()][(idx+1)%3];
    }
    int v_tgt(std::vector<Vector3i> &faces){
        return faces[face()][(idx+2)%3];
    }
};

struct HalfedgeMesh {
    std::vector<Vector3d> V;
    std::vector<Vector3i> F;
    std::vector<Halfedge> hEList;
    std::vector<int> h_out;

    void read_file(std::string& filename) {
        std::ifstream ifs(filename);
        if (ifs.fail()){
            std::cerr << "Failed to open file." << "\n";
            std::exit(1);
        }
        std::string line;
        if (std::filesystem::path(filename).extension() == ".obj") {
            while (std::getline(ifs, line)){
                if (line.empty() != 0 || line[0] == '#') {
                    // Do nothing.
                } else if (line[0] == 'v') {
                    Vector3d v;
                    std::istringstream string_in{line.substr(1)};
                    string_in >> v(0) >> v(1) >> v(2);
                    V.push_back(v);
                } else if (line[0] == 'f') {
                    Vector3i f;
                    std::istringstream string_in{line.substr(1)};
                    string_in >> f(0) >> f(1) >> f(2);
                    f = f - Vector3i{1, 1, 1};
                    F.push_back(f);
                }
            }
        } else if (std::filesystem::path(filename).extension() == ".off") {
            std::getline(ifs, line);    // skip
            std::getline(ifs, line);
            int vsize, fsize;
            sscanf(line.data(), "%d %d", &vsize, &fsize);

            for (int i=0; i<vsize; i++) {
                std::getline(ifs, line);
                Vector3d v;
                std::istringstream string_in{line.substr(1)};
                string_in >> v(0) >> v(1) >> v(2);
                V.push_back(v);
            }
            for (int i=0; i<fsize; i++) {
                std::getline(ifs, line);
                Vector3i f;
                std::istringstream string_in{line.substr(1)};
                string_in >> f(0) >> f(1) >> f(2);
                F.push_back(f);
            }
        }
    }

    void make_halfedge_list() {
        std::pair<int, int> key;
        std::pair<int, int> keyswap;
        std::map<std::pair<int, int>, int> map;
        for (int i=0; i<F.size(); i++) {
            for (int j=0; j<3; j++) {
                Halfedge h;
                h.idx = 3*i+j;

                key = std::make_pair(h.v_src(F), h.v_tgt(F));
                keyswap = std::make_pair(h.v_tgt(F),h.v_src(F));
                if (map.contains(keyswap)) {
                    h.h_opp = map.at(keyswap);
                    hEList[map.at(keyswap)].h_opp = 3*i+j;
                } else {
                    h.h_opp = -1;
                    map.emplace(key, 3*i+j);
                }

                hEList.push_back(h);
            }
        }
        // h_out を計算
        h_out.resize(V.size());
        for (int i=0; i<hEList.size(); i++) {
            //  h_out に境界半辺が保存されていない場合のみ更新
            if (hEList[ h_out[hEList[i].v_src(F)] ].h_opp != -1) {
                h_out[hEList[i].v_src(F)] = i;
            }
        }
        // 境界半辺の h_opp に次の境界半辺を保存する
        /*for (int i=0; i<h_out.size(); i++) {
            if(hEList[h_out[i]].h_opp == -1){
                hEList[h_out[i]].h_opp = -h_out[ hEList[h_out[i]].v_tgt(F) ] - 1;
                // h_i1.h_opp = -i2-1;
            }
        }*/
    }
    int h_ccw(int i) {
        if(i < 0) return i;
        return hEList[ hEList[i].h_prev() ].h_opp;
    }
    int h_cw(int i) {
        if (i < 0) return i;
        else if (hEList[i].h_opp < 0) return hEList[i].h_opp;
        return hEList[ hEList[i].h_opp ].h_next();
    }
};
struct Ray {
    Vector3d origin;
    Vector3d direction;
};

struct Sphere {
    Vector3d center;
    double radius;
};

void update_minmax(double v, int id, double& min, double& max, int& min_id, int& max_id) {
    if (v < min) {
        min = v;
        min_id = id;
    }
    if (v > max) {
        max = v;
        max_id = id;
    }
}
Sphere update_minimum_enclosing_sphere(Vector3d c_old, double r_old, Vector3d pi) {
    Sphere s;
    Vector3d c_new;
    double r_new;
    double d = (pi - c_old).norm() - r_old;
    if ( (pi - c_old).norm() > r_old ) {
        c_new = c_old + (d / 2.0) * ( (pi - c_old) / (pi - c_old).norm() );
        r_new = r_old + d / 2.0;
    } else {
        c_new = c_old;
        r_new = r_old;
    }
    s = {c_new, r_new};
    return s;
}
Sphere minimum_enclosing_sphere(HalfedgeMesh& mesh) {
    Sphere s;
    double minx = mesh.V[0](0), miny = mesh.V[0](1), minz = mesh.V[0](2);
    double maxx = mesh.V[0](0), maxy = mesh.V[0](1), maxz = mesh.V[0](2);
    int minx_id=0, miny_id=0, minz_id=0;
    int maxx_id=0, maxy_id=0, maxz_id=0;
    for (int i=0; i<mesh.V.size(); i++) {
        update_minmax(mesh.V[i](0), i, minx, maxx, minx_id, maxx_id);               // (1, 0, 0)方向の最小値，最大値を更新
        update_minmax(mesh.V[i](1), i, miny, maxy, miny_id, maxy_id);               // (0, 1, 0)方向の最小値，最大値を更新
        update_minmax(mesh.V[i](2), i, minz, maxz, minz_id, maxz_id);               // (0, 0, 1)方向の最小値，最大値を更新
    }

    s.center = ( mesh.V[minx_id] + mesh.V[maxx_id] ) / 2.0;                         // (1, 0, 0)方向の最小点，最大点を初期球とする
    s.radius = ( mesh.V[minx_id] - mesh.V[maxx_id] ).norm() / 2.0;
    s = update_minimum_enclosing_sphere(s.center, s.radius, mesh.V[miny_id]);        // (0, 1, 0)方向の最小値を加える
    s = update_minimum_enclosing_sphere(s.center, s.radius, mesh.V[maxy_id]);        // (0, 1, 0)方向の最大値を加える
    s = update_minimum_enclosing_sphere(s.center, s.radius, mesh.V[minz_id]);        // (0, 0, 1)方向の最小値を加える
    s = update_minimum_enclosing_sphere(s.center, s.radius, mesh.V[maxz_id]);        // (0, 0, 1)方向の最大値を加える

    for (int i=0; i<mesh.V.size(); i++) {
        s = update_minimum_enclosing_sphere(s.center, s.radius, mesh.V[i]);
    }
    return s;
}

bool ray_mesh_intersection(HalfedgeMesh& mesh,
                         Ray& r, int& t_idx, Vector3d& bary_coords,
                         std::vector<Vector3d>& normalF) {

    double t_min;
    int flag = 0;   // 最近点のパラメータ t_min を初期化するためのフラグ
    for (int i = 0; i < mesh.F.size(); i++) {
        Vector3i f = mesh.F[i];
        double t = normalF[i].dot(mesh.V[f(0)] - r.origin) / normalF[i].dot(r.direction);     // ray parameter
        if (t >= 0) {   // レイのパラメータ t が非負のとき，平面と交差する
            Vector3d bc_temp;
            Vector3d bc_dash;
            for (int j = 0; j < 3; j++) {
                Matrix3d mat;
                mat.col(0) = r.direction;
                mat.col(1) = mesh.V[f((j+1)%3)] - r.origin;
                mat.col(2) = mesh.V[f((j+2)%3)] - r.origin;
                bc_dash(j) = mat.determinant();
            }
            bc_temp = bc_dash / (bc_dash(0) + bc_dash(1) + bc_dash(2));
            if (bc_temp(0) >= 0 && bc_temp(1) >= 0 && bc_temp(2) >= 0) {    // 重心座標が全て非負のとき，交点は三角形の内部にある
                if (flag == 0) {    // 交差する三角形を初期化
                    t_min = t;
                    t_idx = i;
                    bary_coords = bc_temp;
                    flag = 1;
                }
                if (t < t_min) {    // パラメータ t が最小のとき可視点
                    t_min = t;
                    t_idx = i;
                    bary_coords = bc_temp;
                }
            }
        }
    }
    return t_idx == -1 ? false : true;
}

void compute_RGB(HalfedgeMesh& mesh, Ray& ray,
                 int& r, int& g, int& b,
                 double& ambient, double& diffuse, Vector3d& light,
                 std::vector<Vector3d>& normalV, std::vector<Vector3d>& normalF) {
    r = 255;
    g = 255;
    b = 255;
    int t_idx = -1;     // 光線が最初に交差する面の添え字
    Vector3d bary_coords;
    bool hit = ray_mesh_intersection(mesh, ray, t_idx, bary_coords, normalF);
    if (hit) {
        Vector3d normal = {0.0, 0.0, 0.0};
        normal += bary_coords(0) * normalV[mesh.F[t_idx][0]];
        normal += bary_coords(1) * normalV[mesh.F[t_idx][1]];
        normal += bary_coords(2) * normalV[mesh.F[t_idx][2]];
        normal.normalize();

        double l = ambient + diffuse * std::max(normal.dot(light), 0.0);
        int c = l * 256;
        if (c > 255) c = 255;
        r = c;
        g = c;
        b = c;
    }
}

void draw_mesh(HalfedgeMesh& mesh, int& nx, int& ny,
               std::vector<Vector3d>& normalV, std::vector<Vector3d>& normalF) {
    double ambient = 0.1, diffuse = 0.9;
    Vector3d light = {1.0, 1.0, 1.0};
    light.normalize();

    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < 2*nx; ++ix) {
            double x = -1.0 + (1.0 / nx) * (ix + 0.5);
            double y = -1.0 + (2.0 / ny) * ((ny - 1 - iy) + 0.75);
            double z = 2.0;
            Ray ray;
            ray.origin = {x, y, z};
            ray.direction = {0.0, 0.0, -1.0};

            int r, g, b;
            // 文字色を設定
            compute_RGB(mesh, ray, r, g, b, ambient, diffuse, light, normalV, normalF);
            std::cout << "\033[38;2;" << r << ";" << g << ";" << b;

            // 背景色を設定
            y = -1.0 + (2.0 / ny) * ((ny - 1 - iy) + 0.25);
            ray.origin = {x, y, z};
            compute_RGB(mesh, ray, r, g, b, ambient, diffuse, light, normalV, normalF);
            std::cout << ";48;2;" << r << ";" << g << ";" << b << "m";

            // "▀"　を出力
            std::cout << "\u2580";
        }
        std::cout << "\033[0m" << std::endl;
    }
}


// =====================================================================
// main
// =====================================================================
int main(int argc, char *argv[]){
    auto start = std::chrono::system_clock::now();

    //    メッシュファイル読み込み
    std::string filename;
    if (argc != 4) {
        std::cout << "usage : " << argv[0] << "  filename  width  height" << std::endl;
        std::exit(1);
    }
    filename = std::string(argv[1]);

    HalfedgeMesh mesh;
    mesh.read_file(filename);
    mesh.make_halfedge_list();

    // メッシュを正規化
    Sphere s = minimum_enclosing_sphere(mesh);    // 最小包含球の計算
    for (auto& v : mesh.V) {
        v = (v - s.center) / s.radius;
    }
    // 頂点法線・面法線を計算
    std::vector<Vector3d> normalV(mesh.V.size(), Vector3d::Zero());
    std::vector<Vector3d> normalF(mesh.F.size(), Vector3d::Zero());
    for (int i = 0; i < mesh.F.size(); i++) {
        Vector3i f = mesh.F[i];
        normalF[i] = ( (mesh.V[f(1)] - mesh.V[f(0)]).cross(mesh.V[f(2)] - mesh.V[f(0)]) ).normalized();
        normalV[f[0]] += normalF[i];
        normalV[f[1]] += normalF[i];
        normalV[f[2]] += normalF[i];
    }
    for (auto& vn : normalV) {
        vn.normalize();
    }

    // 描画
    int nx = std::stoi(std::string(argv[2]));
    int ny = std::stoi(std::string(argv[3]));
    draw_mesh(mesh, nx, ny, normalV, normalF);

    // 出力
    auto end = std::chrono::system_clock::now();
    using namespace std::chrono_literals;
    std::cerr << "実行時間 : " << (end - start) / 1.0s << " 秒\n";
}
