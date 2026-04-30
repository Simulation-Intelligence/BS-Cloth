#include "mesh.h"
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <set>

namespace IPC
{

unsigned long long calculate_Triangle_hash(const uint64_t &index0, const uint64_t &index1, const uint64_t &index2,
                                           const uint64_t &length) {
    unsigned long long hashVal = (index2 * length + index1) * length + index0;
    return hashVal;
}

class Triangle {

public:
    uint64_t key[3];

    Triangle(const uint64_t *p_key) {
        key[0] = p_key[0];
        key[1] = p_key[1];
        key[2] = p_key[2];
    }

    Triangle(uint64_t key0, uint64_t key1, uint64_t key2) {
        key[0] = key0;
        key[1] = key1;
        key[2] = key2;
    }

    uint64_t operator[](int i) const {
        assert(0 <= i && i <= 2);
        return key[i];
    }

    bool operator<(const Triangle &right) const {
        if (key[0] < right.key[0]) {
            return true;
        } else if (key[0] == right.key[0]) {
            if (key[1] < right.key[1]) {
                return true;
            } else if (key[1] == right.key[1]) {
                if (key[2] < right.key[2]) {
                    return true;
                }
            }
        }
        return false;
    }

    bool operator==(const Triangle &right) const {
        return key[0] == right[0] && key[1] == right[1] && key[2] == right[2];
    }
};

void split(string str, vector<string> &v, string spacer) {
    int pos1, pos2;
    int len = spacer.length();
    pos1 = 0;
    pos2 = str.find(spacer);
    while (pos2 != string::npos) {
        v.push_back(str.substr(pos1, pos2 - pos1));
        pos1 = pos2 + len;
        pos2 = str.find(spacer, pos1);
    }
    if (pos1 != str.length())
        v.push_back(str.substr(pos1));
}

void mesh2D::InitMesh(int type, double scale) {
    if (type == 0) {
        mesh_circle circle;
        for (int i = 0; i < circle.vertexNum; i++) {
            Vector2d vertex = scale * Vector2d(circle.vertex[i][0], circle.vertex[i][1]);
            Vector2d force = Vector2d(0, 0);
            Vector2d velocity = Vector2d(0, 0);
            double mass = 0;
            vertexes.push_back(vertex);
            forces.push_back(force);
            velocities.push_back(velocity);
            masses.push_back(mass);
        }

        for (int i = 0; i < circle.triangleNum; i++) {
            Vector3i triangle;
            triangle[0] = (circle.index[i][0]);
            triangle[1] = (circle.index[i][1]);
            triangle[2] = (circle.index[i][2]);
            //Vector3i triangle = Vector3i(circle.index[i][0], circle.index[i][1], circle.index[i][2]);
            triangles.push_back(triangle);
        }

        triangleNum = circle.triangleNum;
        vertexNum = circle.vertexNum;
    } else if (type == 1) {
        mesh_rectangle rectangle;
        for (int i = 0; i < rectangle.vertexNum; i++) {
            Vector2d vertex = scale * Vector2d(rectangle.vertex[i][0], rectangle.vertex[i][1]);
            Vector2d force = Vector2d(0, 0);
            Vector2d velocity = Vector2d(0, 0);
            double mass = 0;
            vertexes.push_back(vertex);
            forces.push_back(force);
            velocities.push_back(velocity);
            masses.push_back(mass);
        }

        for (int i = 0; i < rectangle.triangleNum; i++) {
            Vector3i triangle;
            triangle[0] = (rectangle.index[i][0]);
            triangle[1] = (rectangle.index[i][1]);
            triangle[2] = (rectangle.index[i][2]);
            //Vector3i triangle = Vector3i(circle.index[i][0], circle.index[i][1], circle.index[i][2]);
            triangles.push_back(triangle);
        }

        triangleNum = rectangle.triangleNum;
        vertexNum = rectangle.vertexNum;
    }
}

void mesh3D::InitMesh(int type, double scale) {
    mesh_cuboid cuboid;
    for (int i = 0; i < cuboid.vertexNum; i++) {
        Vector3d vertex = scale * Vector3d(cuboid.vertex[i][0], cuboid.vertex[i][1], cuboid.vertex[i][2]);
        Vector3d force = Vector3d(0, 0, 0);
        Vector3d velocity = Vector3d(0, 0, 0);
        Vector3d d_velocity = Vector3d(0, 0, 0);
        Vector3d d_pos = Vector3d(0, 0, 0);
        double mass = 0;
        Matrix3d Constraint;
        Constraint.setIdentity();
        vertices.push_back(vertex);
        //forces.push_back(force);
        velocities.push_back(velocity);
        //d_velocities.push_back(d_velocity);
        Constraints.push_back(Constraint);
        masses.push_back(mass);
        //d_positions.push_back(d_pos);
    }

    for (int i = 0; i < cuboid.tetrahedralNum; i++) {
        Vector4i tetrahedra;
        for (int j = 0; j < 4; j++) {
            tetrahedra[j] = (cuboid.index[i][0]);
        }
        tetrahedrals.push_back(tetrahedra);
    }

    tetrahedralNum = cuboid.tetrahedralNum;
    vertexNum = cuboid.vertexNum;
}

void mesh3D::load_test(double scale, int num) {
    vertexNum = 8 * num;
    tetrahedralNum = 5 * num;
    mesh_cuboid test;
    double xmin = 1e32, ymin = 1e32, zmin = 1e32;
    double xmax = -1e32, ymax = -1e32, zmax = -1e32;

    for (int i = 0; i < 5; i++) {
        Vector4i tet;
        for (int j = 0; j < 4; j++) {
            tet[j] = (test.index[i][j] + 8);
        }
        tetrahedrals.push_back(tet);

        //tetra_fiberDir.push_back(Vector3d(0, 1, 0));
        //rehabilitate.push_back(-1);
        //isInside.push_back(-1);
        //isJaw.push_back(false);
        //isSkull.push_back(false);
        //isMouth.push_back(false);
        //isMuscle.push_back(false);
    }

    for (int i = 0; i < 8; i++) {
        Vector3d vert =
                scale * Vector3d(test.vertex[i][0], test.vertex[i][1], test.vertex[i][2]) + Vector3d(0.0, 0.0, 0);
        vertices.push_back(vert);

        Vector3d d_velocity = Vector3d(0, 0, 0);
        Matrix3d Constraint;
        Constraint.setIdentity();
        Vector3d force = Vector3d(0, 0, 0);
        Vector3d velocity = Vector3d(0, 0, 0);
        Vector3d d_pos = Vector3d(0, 0, 0);
        double mass = 0;
        //forces.push_back(force);
        velocities.push_back(velocity);
        Constraints.push_back(Constraint);
        //d_velocities.push_back(d_velocity);
        masses.push_back(mass);
        //isDelete.push_back(false);
        //d_positions.push_back(d_pos);
        //externalForce.push_back(Vector3d(0, 0, 0));

        Vector3d pos = vert;
        if (xmin > pos[0]) xmin = pos[0];
        if (ymin > pos[1]) ymin = pos[1];
        if (zmin > pos[2]) zmin = pos[2];
        if (xmax < pos[0]) xmax = pos[0];
        if (ymax < pos[1]) ymax = pos[1];
        if (zmax < pos[2]) zmax = pos[2];
    }
    if (num == 2) {
        for (int i = 0; i < 8; i++) {
            Vector3d vert =
                    scale * Vector3d(test.vertex[i][0], test.vertex[i][1], test.vertex[i][2]) - Vector3d(0.0, 0.5, 0);
            vertices.push_back(vert);

            Vector3d d_velocity = Vector3d(0, 0, 0);
            Matrix3d Constraint;
            Constraint.setIdentity();
            Vector3d force = Vector3d(0, 0, 0);
            Vector3d velocity = Vector3d(0, 0, 0);
            Vector3d d_pos = Vector3d(0, 0, 0);
            double mass = 0;
            //forces.push_back(force);
            velocities.push_back(velocity);
            Constraints.push_back(Constraint);
            //d_velocities.push_back(d_velocity);
            masses.push_back(mass);
            //isDelete.push_back(false);
            //d_positions.push_back(d_pos);
            //externalForce.push_back(Vector3d(0, 0, 0));

            Vector3d pos = vert;
            if (xmin > pos[0]) xmin = pos[0];
            if (ymin > pos[1]) ymin = pos[1];
            if (zmin > pos[2]) zmin = pos[2];
            if (xmax < pos[0]) xmax = pos[0];
            if (ymax < pos[1]) ymax = pos[1];
            if (zmax < pos[2]) zmax = pos[2];
        }

        for (int i = 0; i < 5; i++) {
            Vector4i tet;
            for (int j = 0; j < 4; j++) {
                tet[j] = (test.index[i][j] + 8);
            }
            tetrahedrals.push_back(tet);

            //tetra_fiberDir.push_back(Vector3d(0, 1, 0));
            //rehabilitate.push_back(-1);
            //isInside.push_back(-1);
            //isJaw.push_back(false);
            //isSkull.push_back(false);
            //isMouth.push_back(false);
            //isMuscle.push_back(false);
        }
    }
    minCorner = Vector3d(xmin, ymin, zmin);
    maxCorner = Vector3d(xmax, ymax, zmax);
    V_prev = vertices;
}

bool mesh3D::load_tetrahedraMesh(const ::std::string &filename, double scale, Vector3d position_offset) {

    ifstream ifs(filename);
    if (!ifs) {

        fprintf(stderr, "unable to read file %s\n", filename.c_str());
        ifs.close();
        exit(-1);
        return false;
    }

    double x, y, z;
    int index0, index1, index2, index3;
    string line = "";
    int nodeNumber = 0;
    int elementNumber = 0;
    int offset = vertices.size();
    Vector3d tempMinConer, tempMaxConer;
    // int cnt = 0;
    while (getline(ifs, line)) {
        if (line.length() <= 1) continue;
        // if (cnt < 20) {
        // 	printf("line [%d]: (%d) [%s]\n", cnt, (int)line.length(), line.data());
        // }
        if (line.substr(1, 5) == "Nodes") {
            getline(ifs, line);
            nodeNumber = atoi(line.c_str());
            vertexNum = nodeNumber;

            // printf("%d nodes in total\n", vertexNum);
            double xmin = 1e32, ymin = 1e32, zmin = 1e32;
            double xmax = -1e32, ymax = -1e32, zmax = -1e32;
            for (int i = 0; i < nodeNumber; i++) {
                getline(ifs, line);
                vector<::std::string> nodePos;
                ::std::string spacer = " ";
                split(line, nodePos, spacer);
                x = atof(nodePos[1].c_str());
                y = atof(nodePos[2].c_str());
                z = atof(nodePos[3].c_str());
                Vector3d d_velocity = Vector3d(0, 0, 0);
                Vector3d vertex = scale * Vector3d(x, y, z) + position_offset;
                Matrix3d Constraint;
                Constraint.setIdentity();
                //Vector3d force = Vector3d(0, 0, 0);
                Vector3d velocity = Vector3d(0, 0, 0);
                //Vector3d d_pos = Vector3d(0, 0, 0);
                double mass = 0;
                vertices.push_back(vertex);
                //forces.push_back(force);
                velocities.push_back(velocity);
                Constraints.push_back(Constraint);
                //d_velocities.push_back(d_velocity);
                masses.push_back(mass);
                //isDelete.push_back(false);
                //d_positions.push_back(d_pos);
                //externalForce.push_back(Vector3d(0, 0, 0));
                int boundaryType = 0;
                boundaryTypes.push_back(boundaryType);
                Vector3d pos = vertex;
                if (xmin > pos[0]) xmin = pos[0];
                if (ymin > pos[1]) ymin = pos[1];
                if (zmin > pos[2]) zmin = pos[2];
                if (xmax < pos[0]) xmax = pos[0];
                if (ymax < pos[1]) ymax = pos[1];
                if (zmax < pos[2]) zmax = pos[2];
            }
            tempMinConer = Vector3d(xmin, ymin, zmin);
            tempMaxConer = Vector3d(xmax, ymax, zmax);

            if (maxCorner[0] < tempMaxConer[0]) maxCorner[0] = tempMaxConer[0];
            if (maxCorner[1] < tempMaxConer[1]) maxCorner[1] = tempMaxConer[1];
            if (maxCorner[2] < tempMaxConer[2]) maxCorner[2] = tempMaxConer[2];
            if (minCorner[0] > tempMinConer[0]) minCorner[0] = tempMinConer[0];
            if (minCorner[1] > tempMinConer[1]) minCorner[1] = tempMinConer[1];
            if (minCorner[2] > tempMinConer[2]) minCorner[2] = tempMinConer[2];

        }

        if (line.substr(1, 8) == "Elements") {
            getline(ifs, line);
            elementNumber = atoi(line.c_str());
            tetrahedralNum = elementNumber;
            for (int i = 0; i < elementNumber; i++) {
                getline(ifs, line);

                vector<::std::string> elementIndexex;
                ::std::string spacer = " ";
                split(line, elementIndexex, spacer);
                index0 = atoi(elementIndexex[3].c_str()) - 1;
                index1 = atoi(elementIndexex[4].c_str()) - 1;
                index2 = atoi(elementIndexex[5].c_str()) - 1;
                index3 = atoi(elementIndexex[6].c_str()) - 1;

                Vector4i tetrahedra;
                tetrahedra[0] = (index0 + offset);
                tetrahedra[1] = (index1 + offset);
                tetrahedra[2] = (index2 + offset);
                tetrahedra[3] = (index3 + offset);

                tetrahedrals.push_back(tetrahedra);

                //tetra_fiberDir.push_back(Vector3d(0, 1, 0));
                //rehabilitate.push_back(-1);
                //isInside.push_back(-1);
                //isJaw.push_back(false);
                //isSkull.push_back(false);
                //isMouth.push_back(false);
                //isMuscle.push_back(false);
            }
            break;
        }
        // cnt++;
    }

    if ((tempMaxConer - tempMinConer).norm() > (objMaxConer - objMinConer).norm()) {
        objMaxConer = tempMaxConer;
        objMinConer = tempMinConer;
    }

    ifs.close();
    V_prev = vertices;

    D12x12Num = 0;
    D9x9Num = 0;
    D6x6Num = 0;
    D3x3Num = 0;


    return true;
}


bool mesh3D::load_tetrahedraMesh_IPC_TetMesh(const ::std::string &filename, double scale, Vector3d position_offset) {
    ifstream ifs(filename);
    if (!ifs) {

        fprintf(stderr, "unable to read file %s\n", filename.c_str());
        ifs.close();
        exit(-1);
        return false;
    }

    double x, y, z;
    int index0, index1, index2, index3;
    string line = "";
    int nodeNumber = 0;
    int elementNumber = 0;
    int offset = vertices.size();
    printf("offset %d\n", offset);
    Vector3d tempMinConer, tempMaxConer;
    while (getline(ifs, line)) {
        if (line.length() <= 1) continue;
        if (line == "$Nodes") {
            getline(ifs, line);
            nodeNumber = atoi(line.c_str());

            double xmin = 1e32, ymin = 1e32, zmin = 1e32;
            double xmax = -1e32, ymax = -1e32, zmax = -1e32;
            for (int i = 0; i < nodeNumber; i++) {
                getline(ifs, line);
                vector<::std::string> nodePos;
                ::std::string spacer = " ";
                split(line, nodePos, spacer);
                x = atof(nodePos[0].c_str());
                y = atof(nodePos[1].c_str());
                z = atof(nodePos[2].c_str());
                Vector3d d_velocity = Vector3d(0, 0, 0);
                Vector3d vertex = Vector3d(scale * x - position_offset.x(), scale * y - position_offset.y(),
                                           scale * z - position_offset.z());
                Matrix3d Constraint;
                Constraint.setIdentity();
                Vector3d force = Vector3d(0, 0, 0);
                Vector3d velocity = Vector3d(0, 0, 0);
                Vector3d d_pos = Vector3d(0, 0, 0);
                double mass = 0;
                vertices.push_back(vertex);
                //forces.push_back(force);
                velocities.push_back(velocity);
                Constraints.push_back(Constraint);
                //d_velocities.push_back(d_velocity);
                masses.push_back(mass);
                //isDelete.push_back(false);
                //d_positions.push_back(d_pos);
                //externalForce.push_back(Vector3d(0, 0, 0));

                int boundaryType = 0;
                boundaryTypes.push_back(boundaryType);
                Vector3d pos = vertex;
                if (xmin > pos[0]) xmin = pos[0];
                if (ymin > pos[1]) ymin = pos[1];
                if (zmin > pos[2]) zmin = pos[2];
                if (xmax < pos[0]) xmax = pos[0];
                if (ymax < pos[1]) ymax = pos[1];
                if (zmax < pos[2]) zmax = pos[2];
            }
            tempMinConer = Vector3d(xmin, ymin, zmin);
            tempMaxConer = Vector3d(xmax, ymax, zmax);

            if (maxCorner[0] < tempMaxConer[0]) maxCorner[0] = tempMaxConer[0];
            if (maxCorner[1] < tempMaxConer[1]) maxCorner[1] = tempMaxConer[1];
            if (maxCorner[2] < tempMaxConer[2]) maxCorner[2] = tempMaxConer[2];
            if (minCorner[0] > tempMinConer[0]) minCorner[0] = tempMinConer[0];
            if (minCorner[1] > tempMinConer[1]) minCorner[1] = tempMinConer[1];
            if (minCorner[2] > tempMinConer[2]) minCorner[2] = tempMinConer[2];
        }

        if (line == "$Elements") {
            getline(ifs, line);
            elementNumber = atoi(line.c_str());
            for (int i = 0; i < elementNumber; i++) {
                getline(ifs, line);

                vector<::std::string> elementIndexex;
                ::std::string spacer = " ";
                split(line, elementIndexex, spacer);
                index0 = atoi(elementIndexex[1].c_str()) - 1;
                index1 = atoi(elementIndexex[2].c_str()) - 1;
                index2 = atoi(elementIndexex[3].c_str()) - 1;
                index3 = atoi(elementIndexex[4].c_str()) - 1;

                //vector<uint64_t> tetrahedra;
                //tetrahedra.push_back(index0 + offset);
                //tetrahedra.push_back(index1 + offset);
                //tetrahedra.push_back(index2 + offset);
                //tetrahedra.push_back(index3 + offset);
                //tetrahedras.push_back(tetrahedra);

                Vector4i tetrahedra;
                tetrahedra[0] = (index0 + offset);
                tetrahedra[1] = (index1 + offset);
                tetrahedra[2] = (index2 + offset);
                tetrahedra[3] = (index3 + offset);

                tetrahedrals.push_back(tetrahedra);
                //tetra_fiberDir.push_back(Vector3d(0, 1, 0));
                //rehabilitate.push_back(-1);
                //isInside.push_back(-1);
                //isJaw.push_back(false);
                //isSkull.push_back(false);
                //isMouth.push_back(false);
                //isMuscle.push_back(false);
            }
            break;
        }
    }
    ifs.close();
    V_prev = vertices;
    vertexNum = vertices.size();
    tetrahedralNum = tetrahedrals.size();

    if ((tempMaxConer - tempMinConer).norm() > (objMaxConer - objMinConer).norm()) {
        objMaxConer = tempMaxConer;
        objMinConer = tempMinConer;
    }

    D12x12Num = 0;
    D9x9Num = 0;
    D6x6Num = 0;
    D3x3Num = 0;
    return true;

}

void mesh3D::getSurface() {
    surfBackLinks.resize(this->vertices.size());
    
    // detecting surfaces of tetrahedra
   
    //uint64_t length = vertexNum;
    //auto triangle_hash = [&](const Triangle &tri) {
    //    return length * (length * tri[0] + tri[1]) + tri[2];
    //};
    ////vector<Vector4i> surface;
    //::std::unordered_map<Triangle, uint64_t, decltype(triangle_hash)> tri2Tet(4 * tetrahedralNum, triangle_hash);
    //for (int i = 0; i < tetrahedralNum; i++) {
    //    const auto &triI = tetrahedrals[i]; 
    //    for (int j = 0; j < 4; j++) {
    //        const Triangle triVInd = Triangle(triI[j % 4], triI[(1 + j) % 4], triI[(2 + j) % 4]);
    //        if (tri2Tet.find(Triangle(triVInd[0], triVInd[1], triVInd[2])) != tri2Tet.end()) {
    //            tri2Tet[Triangle(triVInd[0], triVInd[1], triVInd[2])] = tetrahedralNum + 1;
    //        } else if (tri2Tet.find(Triangle(triVInd[0], triVInd[2], triVInd[1])) != tri2Tet.end()) {
    //            tri2Tet[Triangle(triVInd[0], triVInd[2], triVInd[1])] = tetrahedralNum + 1;
    //        } else if (tri2Tet.find(Triangle(triVInd[1], triVInd[0], triVInd[2])) != tri2Tet.end()) {
    //            tri2Tet[Triangle(triVInd[1], triVInd[0], triVInd[2])] = tetrahedralNum + 1;
    //        } else if (tri2Tet.find(Triangle(triVInd[1], triVInd[2], triVInd[0])) != tri2Tet.end()) {
    //            tri2Tet[Triangle(triVInd[1], triVInd[2], triVInd[0])] = tetrahedralNum + 1;
    //        } else if (tri2Tet.find(Triangle(triVInd[2], triVInd[0], triVInd[1])) != tri2Tet.end()) {
    //            tri2Tet[Triangle(triVInd[2], triVInd[0], triVInd[1])] = tetrahedralNum + 1;
    //        } else if (tri2Tet.find(Triangle(triVInd[2], triVInd[1], triVInd[0])) != tri2Tet.end()) {
    //            tri2Tet[Triangle(triVInd[2], triVInd[1], triVInd[0])] = tetrahedralNum + 1;
    //        } else {
    //            tri2Tet[Triangle(triVInd[0], triVInd[1], triVInd[2])] = i;
    //        }
    //    }
    //}

    //for (const auto &triI : tri2Tet) {
    //    const uint64_t &tetId = triI.second;
    //    const Triangle &triVInd = triI.first;
    //    if (tetId < tetrahedralNum) {
    //        Vector3d vec1 = vertices[triVInd[1]] - vertices[triVInd[0]];
    //        Vector3d vec2 = vertices[triVInd[2]] - vertices[triVInd[0]];
    //        int id3 = 0;
    //        for (int i = 0; i < 4; i++) {
    //            if (tetrahedrals[tetId][i] != triVInd[0]
    //                && tetrahedrals[tetId][i] != triVInd[1]
    //                && tetrahedrals[tetId][i] != triVInd[2]) {
    //                id3 = tetrahedrals[tetId][i];
    //                break;
    //            } 
    //        }

    //        Vector3d vec3 = vertices[id3] - vertices[triVInd[0]];
    //        Vector3d n = vec1.cross(vec2);
    //        if (n.dot(vec3) < 0) {
    //            size_t curSurfaceIdx = surface.size();
    //            surfBackLinks[triVInd[0]].surfaces.push_back(curSurfaceIdx);
    //            surfBackLinks[triVInd[1]].surfaces.push_back(curSurfaceIdx);
    //            surfBackLinks[triVInd[2]].surfaces.push_back(curSurfaceIdx);
    //            surface.push_back(Vector4i(triVInd[0], triVInd[1], triVInd[2], tetId));
    //        } else {
    //            size_t curSurfaceIdx = surface.size();
    //            surfBackLinks[triVInd[0]].surfaces.push_back(curSurfaceIdx);
    //            surfBackLinks[triVInd[1]].surfaces.push_back(curSurfaceIdx);
    //            surfBackLinks[triVInd[2]].surfaces.push_back(curSurfaceIdx);
    //            surface.push_back(Vector4i(triVInd[0], triVInd[2], triVInd[1], tetId));
    //        }
    //    }
    //}

    unsigned int patchCnt = this->patchSurfaceSep.size();

    for(const auto & tri : triangles){
        size_t curSurfaceIdx = surface.size();
        surfBackLinks[tri[0]].surfaces.push_back(curSurfaceIdx);
        surfBackLinks[tri[1]].surfaces.push_back(curSurfaceIdx);
        surfBackLinks[tri[2]].surfaces.push_back(curSurfaceIdx);
        surface.emplace_back(Vector4i(tri[0], tri[1], tri[2], 0));
    }
    //vector<bool> flag(vertexNum, false);
    //for (const auto &cTri: surface) {
    //    for (int i = 0; i < 3; i++) {
    //        if (!flag[cTri[i]]) {
    //            size_t curVertIdx = surfVerts.size();
    //            surfBackLinks[cTri[i]].surfVerts.push_back(curVertIdx);
    //            surfVerts.push_back(cTri[i]);
    //            flag[cTri[i]] = true;
    //        }
    //    }
    //}
    for (size_t i = 0; i != this->vertices.size(); ++i) {
        surfVerts.push_back(i);
        surfBackLinks[i].surfVerts.push_back(i);
    }

    set<pair<uint64_t, uint64_t>> SFEdges_set;

    // this->patchEdgeSep[0] = 0;
    // for (unsigned int patchIdx = 0; patchIdx != this->patchSurfaceSep.size() - 1; ++patchIdx) {
    //     unsigned int patchSurfStartIdx = this->patchSurfaceSep[patchIdx];
    //     unsigned int patchSurfEndIdx = (patchIdx == (this->patchSurfaceSep.size() - 2)) ?
    //         this->surface.size() : this->patchSurfaceSep[patchIdx + 1];

    //     printf("patch %d: surface start: %d, end: %d\n", patchIdx, patchSurfStartIdx, patchSurfEndIdx);

    //     for (unsigned int surfIdx = patchSurfStartIdx; surfIdx != patchSurfEndIdx; ++surfIdx) {
    //         const auto& cTri = surface[surfIdx];
    //         for (int i = 0; i < 3; i++) {
    //             if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[1], cTri[0])) == SFEdges_set.end()) {
    //                 SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[0], cTri[1]));
    //             }
    //             if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[2], cTri[1])) == SFEdges_set.end()) {
    //                 SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[1], cTri[2]));
    //             }
    //             if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[0], cTri[2])) == SFEdges_set.end()) {
    //                 SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[2], cTri[0]));
    //             }
    //         }
    //     }

    //     for (const auto& edge : SFEdges_set) {
    //         size_t curEdgeIdx = surfEdges.size();
    //         surfBackLinks[edge.first].surfEdges.push_back(curEdgeIdx);
    //         surfBackLinks[edge.second].surfEdges.push_back(curEdgeIdx);
    //         surfEdges.push_back(edge);
    //     }
    //     SFEdges_set.clear();
    //     this->patchEdgeSep[patchIdx + 1] = surfEdges.size();
    // }
    
    this->patchEdgeSep[0] = 0;
    for (unsigned int patchIdx = 0; patchIdx != this->patchSurfaceSep.size() - 1; ++patchIdx) {
        unsigned int patchSurfStartIdx = this->patchSurfaceSep[patchIdx];
        unsigned int patchSurfEndIdx = this->patchSurfaceSep[patchIdx + 1];

        printf("patch %d: surface start: %d, end: %d\n", patchIdx, patchSurfStartIdx, patchSurfEndIdx);

        for (unsigned int surfIdx = patchSurfStartIdx; surfIdx != patchSurfEndIdx; ++surfIdx) {
            const auto& cTri = surface[surfIdx];
            for (int i = 0; i < 3; i++) {
                if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[1], cTri[0])) == SFEdges_set.end()) {
                    SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[0], cTri[1]));
                }
                if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[2], cTri[1])) == SFEdges_set.end()) {
                    SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[1], cTri[2]));
                }
                if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[0], cTri[2])) == SFEdges_set.end()) {
                    SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[2], cTri[0]));
                }
            }
        }

        for (const auto& edge : SFEdges_set) {
            size_t curEdgeIdx = surfEdges.size();
            surfBackLinks[edge.first].surfEdges.push_back(curEdgeIdx);
            surfBackLinks[edge.second].surfEdges.push_back(curEdgeIdx);
            surfEdges.push_back(edge);
        }
        SFEdges_set.clear();
        this->patchEdgeSep[patchIdx + 1] = surfEdges.size();
    }

    for (unsigned int surfIdx = this->patchSurfaceSep[this->patchSurfaceSep.size() - 1]; surfIdx != this->surface.size(); ++surfIdx) {
        const auto& cTri = surface[surfIdx];
        for (int i = 0; i < 3; i++) {
            if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[1], cTri[0])) == SFEdges_set.end()) {
                SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[0], cTri[1]));
            }
            if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[2], cTri[1])) == SFEdges_set.end()) {
                SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[1], cTri[2]));
            }
            if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[0], cTri[2])) == SFEdges_set.end()) {
                SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[2], cTri[0]));
            }
        }
    }

    for (const auto& edge : SFEdges_set) {
        size_t curEdgeIdx = surfEdges.size();
        surfBackLinks[edge.first].surfEdges.push_back(curEdgeIdx);
        surfBackLinks[edge.second].surfEdges.push_back(curEdgeIdx);
        surfEdges.push_back(edge);
    }
    SFEdges_set.clear();

    for (unsigned int i = 0; i != patchCnt; ++i) {
        cout << "patch " << i << " :: vertStart: " << this->patchVertexSep[i] 
            << " edgeStart: " << this->patchEdgeSep[i] 
            << " surfStart: " << this->patchSurfaceSep[i] << endl;
    }

    //for (const auto &cTri: surface) {
    //    for (int i = 0; i < 3; i++) {
    //        if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[1], cTri[0])) == SFEdges_set.end()) {
    //            SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[0], cTri[1]));
    //        }
    //        if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[2], cTri[1])) == SFEdges_set.end()) {
    //            SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[1], cTri[2]));
    //        }
    //        if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[0], cTri[2])) == SFEdges_set.end()) {
    //            SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[2], cTri[0]));
    //        }
    //    }
    //}

    //for (const auto &edge : SFEdges_set) {
    //    size_t curEdgeIdx = surfEdges.size();
    //    surfBackLinks[edge.first].surfEdges.push_back(curEdgeIdx);
    //    surfBackLinks[edge.second].surfEdges.push_back(curEdgeIdx);
    //    surfEdges.push_back(edge);
    //}

    //cout << "started test";

    //// temporary tests
    //unsigned int vertexCnt = this->vertices.size();
    //for (unsigned int i = 0; i != vertexCnt; ++i)
    //{
    //    const SurfaceInfo& info = surfBackLinks[i];
    //    cout << "## vertex " << i << " : " << endl;

    //    for (unsigned int idx : info.surfVerts)
    //        cout << "   - surf vert : " << idx << " (global id : " << surfVerts[idx] << ")" << endl;

    //    for (unsigned int idx : info.surfEdges)
    //        cout << "   - surf edge : " << idx << " (global id : " << surfEdges[idx].first << ", " << surfEdges[idx].second << ")" << endl;

    //    for (const auto& edge : this->surfEdges)
    //        if (edge.first == i || edge.second == i)
    //            cout << "       * edge : (global id : " << edge.first << ", " << edge.second << ")" << endl;

    //    for (unsigned int idx : info.surfaces)
    //        cout << "   - surf tri  : " << idx << " (global id : " << surface[idx][0] << ", " << surface[idx][1] << ", " << surface[idx][2] << ")" << endl;

    //    for (const auto& tri : this->surface)
    //        if (tri[0] == i || tri[1] == i || tri[2] == i)
    //            cout << "       * tri : (global id : " << tri[0] << ", " << tri[1] << ", " << tri[2] << ")" << endl;
    //}

    //cout << "finished test";
}

bool mesh3D::output_tetrahedraMesh(const ::std::string &filename) {
    ::std::ofstream outmsh1(filename);
    outmsh1 << "$Nodes\n";
    outmsh1 << vertexNum << endl;
    for (int i = 0; i < vertexNum; i++) {
        outmsh1 << i + 1 << " " << vertices[i][0] << " " <<
                vertices[i][1] << " " <<
                vertices[i][2] << endl;
    }
    outmsh1 << "$Elements\n";
    outmsh1 << tetrahedralNum << endl;
    for (int i = 0; i < tetrahedralNum; i++) {
        outmsh1 << i + 1 << " 4 0 " << tetrahedrals[i][0] << " " <<
                tetrahedrals[i][1] << " " <<
                tetrahedrals[i][2] << " " <<
                tetrahedrals[i][3] << endl;
    }
    outmsh1.close();
    return true;
}

bool mesh3D::load_tetTempData() {
    ::std::string fileVertex = "tempData/vertex.txt";
    ::std::string fileXtileVertex = "tempData/vertexXtile.txt";
    ifstream ifs1(fileVertex);
    ifstream ifs2(fileXtileVertex);

    if (!ifs1 || !ifs2) {
        ifs1.close();
        ifs2.close();
        return false;
    }
    double x, y, z;
    int index = 0;
    while (ifs1 >> x >> y >> z) {
        vertices[index] = Vector3d(x, y, z);
        V_prev[index] = vertices[index];
        index++;
    }
    index = 0;
    while (ifs2 >> x >> y >> z) {
        xTilta[index++] = Vector3d(x, y, z);
    }
    ifs1.close();
    ifs2.close();
    return true;
}

bool mesh3D::output_tetTempData() {
    ::std::string fileVertex = "tempData/vertex.txt";
    ::std::string vertexXtile = "tempData/vertexXtile.txt";


    ::std::ofstream outmsh1(fileVertex);
    for (int i = 0; i < vertexNum; i++) {
        outmsh1 << vertices[i][0] << " " <<
                vertices[i][1] << " " <<
                vertices[i][2] << endl;
    }
    outmsh1.close();

    ::std::ofstream outmsh2(vertexXtile);
    for (int i = 0; i < vertexNum; i++) {
        outmsh2 << xTilta[i][0] << " " <<
                xTilta[i][1] << " " <<
                xTilta[i][2] << endl;
    }
    outmsh2.close();
    return true;
}

bool mesh3D::load_triangleMesh(const string &filename, double scale, Vector3d position_offset, int type) {
    int offset = vertices.size();
    ifstream ifs(filename);
    if (!ifs) {

        fprintf(stderr, "unable to read file %s\n", filename.c_str());
        ifs.close();
        exit(-1);
        return false;
    }
    char buffer[1024];
    string line = "";
    int nodeNumber = 0;
    int elementNumber = 0;
    double x, y, z;

    double xmin = 1e32, ymin = 1e32, zmin = 1e32;
    double xmax = -1e32, ymax = -1e32, zmax = -1e32;
    while (getline(ifs, line)) {
        if (line.length() <= 1) continue;
        string key = line.substr(0, 2);
        stringstream ss(line.substr(2));
        if (key == "v ") {
            ss >> x >> y >> z;
            Vector3d d_velocity = Vector3d(0, 0, 0);
            Vector3d vertex = scale * Vector3d(x, y, z) + position_offset;
            Matrix3d Constraint;
            Constraint.setIdentity();
            Vector3d velocity = Vector3d(0, 0, 0);
            double mass = 0;
            int boundaryType = 0;

            //if(type==2){
            //    //mass = 10;
            //    Constraint.setZero();
            //    boundaryType = 2;
            //}

            vertices.push_back(vertex);
            velocities.push_back(velocity);
            Constraints.push_back(Constraint);
            masses.push_back(mass);
            boundaryTypes.push_back(boundaryType);

            Vector3d pos = vertex;
            if (xmin > pos[0]) xmin = pos[0];
            if (ymin > pos[1]) ymin = pos[1];
            if (zmin > pos[2]) zmin = pos[2];
            if (xmax < pos[0]) xmax = pos[0];
            if (ymax < pos[1]) ymax = pos[1];
            if (zmax < pos[2]) zmax = pos[2];

        } 
        else if (key == "vn") {
            ss >> x >> y >> z;
            Vector3d normal = Vector3d(x, y, z);
//            normals.push_back(normal);
        } 
        else if (key == "f ") {
            if (line.length() >= 1024) {
                printf("[WARN]: skip line due to exceed max buffer length (1024).\n");
                continue;
            }

            ::std::vector<string> fs;

            {
                string buf;
                stringstream ss(line);
                vector<string> tokens;
                while (ss >> buf)
                    tokens.push_back(buf);

                for (size_t index = 3; index < tokens.size(); index += 1) {
                    fs.push_back("f " + tokens[1] + " " + tokens[index - 1] + " " + tokens[index]);
                }
            }

            int uv0, uv1, uv2;

            for (const auto &f: fs) {
                memset(buffer, 0, sizeof(char) * 1024);
                std::copy(f.begin(), f.end(), buffer);

                Vector3i faceVertIndex;
                Vector3i faceNormalIndex;

                if (sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d", &faceVertIndex[0], &uv0, &faceNormalIndex[0],
                           &faceVertIndex[1], &uv1, &faceNormalIndex[1],
                           &faceVertIndex[2], &uv2, &faceNormalIndex[2]) == 9) {
                    faceVertIndex[0] -= 1;
                    faceVertIndex[1] -= 1;
                    faceVertIndex[2] -= 1;

                    faceVertIndex[0] += offset;
                    faceVertIndex[1] += offset;
                    faceVertIndex[2] += offset;
                } else if (sscanf(buffer, "f %d %d %d", &faceVertIndex[0], &faceVertIndex[1], &faceVertIndex[2]) == 3) {
                    faceVertIndex[0] -= 1;
                    faceVertIndex[1] -= 1;
                    faceVertIndex[2] -= 1;

                    faceVertIndex[0] += offset;
                    faceVertIndex[1] += offset;
                    faceVertIndex[2] += offset;
                }

                if (type == 3) {
                    surface.push_back(Vector4i(faceVertIndex[0],faceVertIndex[1],faceVertIndex[2],0));
                } else
                    triangles.push_back(faceVertIndex);

            }
        }
    }
    vertexNum = vertices.size();
    Vector3d tempMinConer, tempMaxConer;
    tempMinConer = Vector3d(xmin, ymin, zmin);
    tempMaxConer = Vector3d(xmax, ymax, zmax);
    if (maxCorner[0] < tempMaxConer[0]) maxCorner[0] = tempMaxConer[0];
    if (maxCorner[1] < tempMaxConer[1]) maxCorner[1] = tempMaxConer[1];
    if (maxCorner[2] < tempMaxConer[2]) maxCorner[2] = tempMaxConer[2];
    if (minCorner[0] > tempMinConer[0]) minCorner[0] = tempMinConer[0];
    if (minCorner[1] > tempMinConer[1]) minCorner[1] = tempMinConer[1];
    if (minCorner[2] > tempMinConer[2]) minCorner[2] = tempMinConer[2];









    ::std::set<::std::pair<int, int>> edge_set;
    ::std::map<::std::pair<int, int>, ::std::vector<int>> edge_map;
    ::std::vector<Eigen::Vector2i> my_edges;
    for (auto tri : triangles) {
        auto x = tri[0];
        auto y = tri[1];
        auto z = tri[2];
        if (x < y) {
            edge_set.insert(make_pair(x, y));
            edge_map[make_pair(x, y)].emplace_back(z);
        }
        else {
            edge_set.insert(make_pair(y, x));
            edge_map[make_pair(y, x)].emplace_back(z);
        }

        if (y < z) {
            edge_set.insert(make_pair(y, z));
            edge_map[make_pair(y, z)].emplace_back(x);
        }
        else {
            edge_set.insert(make_pair(z, y));
            edge_map[make_pair(z, y)].emplace_back(x);
        }

        if (x < z) {
            edge_set.insert(make_pair(x, z));
            edge_map[make_pair(x, z)].emplace_back(y);
        }
        else {
            edge_set.insert(make_pair(z, x));
            edge_map[make_pair(z, x)].emplace_back(y);
        }
    }

    ::std::vector<::std::vector<int>> temp_edges_adj_points;
    for (auto p : edge_set) {
        if (edge_map[p].size() != 2)continue;
        my_edges.emplace_back(p.first, p.second);
        temp_edges_adj_points.emplace_back(edge_map[make_pair(p.first, p.second)]);
    }
    for (int i = 0; i < temp_edges_adj_points.size(); i++) {
        tri_edges.emplace_back(Vector2i(my_edges[i].x(), my_edges[i].y()));
        if (temp_edges_adj_points[i].size() == 2)
            tri_edges_adj_points.emplace_back(Vector2i(temp_edges_adj_points[i][0], temp_edges_adj_points[i][1]));
        else
            tri_edges_adj_points.emplace_back(Vector2i(temp_edges_adj_points[i][0], -1));
    }

    return true;
}

bool mesh_obj::load_mesh(const ::std::string &filename, double scale) {
    ifstream ifs(filename);
    if (!ifs) {

        fprintf(stderr, "unable to read file %s\n", filename.c_str());
        ifs.close();
        exit(-1);
        return false;
    }
    char buffer[1024];
    string line = "";
    int nodeNumber = 0;
    int elementNumber = 0;
    double x, y, z;

    while (getline(ifs, line)) {
        string key = line.substr(0, 2);
        stringstream ss(line.substr(2));
        if (key == "v ") {
            ss >> x >> y >> z;
            Vector3d vertex = scale * Vector3d(x, y, z);
            vertexes.push_back(vertex);
        } else if (key == "vn") {
            ss >> x >> y >> z;
            Vector3d normal = Vector3d(x, y, z);
            normals.push_back(normal);
        } else if (key == "f ") {
            if (line.length() >= 1024) {
                printf("[WARN]: skip line due to exceed max buffer length (1024).\n");
                continue;
            }

            ::std::vector<string> fs;

            {
                string buf;
                stringstream ss(line);
                vector<string> tokens;
                while (ss >> buf)
                    tokens.push_back(buf);

                for (size_t index = 3; index < tokens.size(); index += 1) {
                    fs.push_back("f " + tokens[1] + " " + tokens[index - 1] + " " + tokens[index]);
                }
            }

            int uv0, uv1, uv2;

            for (const auto &f: fs) {
                memset(buffer, 0, sizeof(char) * 1024);
                ::std::copy(f.begin(), f.end(), buffer);

                Vector3i faceVertIndex;
                Vector3i faceNormalIndex;

                sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d", &faceVertIndex[0], &uv0, &faceNormalIndex[0],
                       &faceVertIndex[1], &uv1, &faceNormalIndex[1],
                       &faceVertIndex[2], &uv2, &faceNormalIndex[2]);

                faceVertIndex[0] -= 1;
                faceVertIndex[1] -= 1;
                faceVertIndex[2] -= 1;
                faces.push_back(faceVertIndex);
                facenormals.push_back(faceNormalIndex);

            }
        }
    }

    vertexNum = vertexes.size();
    faceNum = faces.size();

    for (int i = 0; i < faceNum; i++) {
        Vector3d A = vertexes[faces[i][1]] - vertexes[faces[i][0]];
        Vector3d B = vertexes[faces[i][2]] - vertexes[faces[i][0]];
        Vector3d normal = (A.cross(B)).normalized();
        normals.push_back(normal);
    }


    return true;
}

bool model_obj::load_model(const ::std::string &filename) {
    ifstream ifs(filename);
    if (!ifs) {

        fprintf(stderr, "unable to read file %s\n", filename.c_str());
        ifs.close();
        exit(-1);
        return false;
    }

    string meshName;
    //int id = 0;
    float scale = 1;
    while (ifs >> meshName) {
        mesh_obj mesh;
        mesh.load_mesh(meshName, scale);
        meshes.push_back(mesh);
        names.push_back(meshName);
    }

    return true;
}

bool fiber_obj::load_model(const ::std::string *filename) {

    ifstream ifsMuscle(filename[0]), ifsTendonIn(filename[1]), ifsTendonOut(filename[2]);
    if (!ifsMuscle || !ifsTendonIn || !ifsTendonOut) {

        fprintf(stderr, "unable to read file %s\n", filename[0].data());
        ifsMuscle.close();
        ifsTendonIn.close();
        ifsTendonOut.close();
        exit(-1);
        return false;
    }

    string meshName;
    float scale = 1;
    while (ifsMuscle >> meshName) {
        mesh_obj mesh;
        mesh.load_mesh(meshName, scale);
        muscles.push_back(mesh);
        //names.push_back(meshName);
    }
    string tendonName;
    while (ifsTendonIn >> tendonName) {
        mesh_obj mesh;
        mesh.load_mesh(tendonName, scale);
        tendonIns.push_back(mesh);
        //names.push_back(meshName);
    }

    while (ifsTendonOut >> tendonName) {
        mesh_obj mesh;
        mesh.load_mesh(tendonName, scale);
        tendonOuts.push_back(mesh);
        //names.push_back(meshName);
    }

    return true;
}

#include<iostream>

void model_tet::calculate_surface() {

    for (auto &tetMesh: mesh3Ds) {
        tetMesh.getSurface();
    }
}

bool model_tet::load_model(const ::std::string &filename, int offset) {


    return true;
}



vector<Eigen::VectorXd> mesh3D::get_dXn1_dXn() {
    return EKF;
}
vector<Eigen::VectorXd> mesh3D::get_dXn1_dVn() {
    vector<Eigen::VectorXd> dXn1_dVn;
    for (int i = 0; i < vertexNum * 3; i++) {
        dXn1_dVn.push_back(IPC_dt * EKF[i]);
    }
    return dXn1_dVn;
}
vector<Eigen::VectorXd> mesh3D::get_dVn1_dXn() {
    vector<Eigen::VectorXd> dVn1_dXn;
    double one_div_dt = 1.0 / IPC_dt;
    for (int i = 0; i < vertexNum * 3; i++) {
        dVn1_dXn.push_back(one_div_dt * EKF[i]);
    }
    return dVn1_dXn;
}
vector<Eigen::VectorXd> mesh3D::get_dVn1_dVn() {
    return EKF;
}

bool AnimFilter(unsigned int patchIdxI, unsigned int patchIdxJ, unsigned int patchCnt)
{
    // box
    bool diffPatchFilter = false;
    //if (patchIdxI >= 0 && patchIdxJ <= 4 && patchIdxI >= 0 && patchIdxJ <= 4)
    //    diffPatchFilter = true;

    // shuffling
    // bool diffPatchFilter = (patchIdxI != patchCnt - 1) && (patchIdxJ != patchCnt - 1) && (patchIdxI != patchIdxJ);
    // if (patchIdxI > patchIdxJ)
    // {
    //     unsigned int tmp = patchIdxI;
    //     patchIdxI = patchIdxJ;
    //     patchIdxJ = tmp;
    // }
    // if (patchIdxI == 0 && patchIdxJ == 1)
    //     diffPatchFilter = false;

    // hiphop
    //if (patchIdxI > patchIdxJ)
    //{
    //    unsigned int tmp = patchIdxI;
    //    patchIdxI = patchIdxJ;
    //    patchIdxJ = tmp;
    //}
    //bool diffPatchFilter = (patchIdxI != patchCnt - 1) && (patchIdxJ != patchCnt - 1);
    ////bool diffPatchFilter = true;
    //if (patchIdxI == 0 || patchIdxI == 1)
    //    if (patchIdxJ == 8 || patchIdxJ == 9 || patchIdxJ == 10 || patchIdxJ == 11)
    //        diffPatchFilter = false;
    //if (patchIdxI == patchIdxJ && patchIdxI != 6 && patchIdxI != 7)
    //    diffPatchFilter = false;

    // dancing
    //if (patchIdxI > patchIdxJ)
    //{
    //    unsigned int tmp = patchIdxI;
    //    patchIdxI = patchIdxJ;
    //    patchIdxJ = tmp;
    //}
    //bool diffPatchFilter = (patchIdxI != patchCnt - 1) && (patchIdxJ != patchCnt - 1) && (patchIdxI != patchIdxJ);
    ////if (patchIdxI == 0 || patchIdxI == 1)
    ////    if (patchIdxJ == 12 || patchIdxJ == 13)
    ////        diffPatchFilter = false;

    //if ((patchIdxI == 8 && patchIdxJ == 10) || (patchIdxI == 9 && patchIdxJ == 11))
    //    diffPatchFilter = false;

    //// In all cases, if the body mesh intersects with itself, ignore the contact:
    if (patchIdxI == patchCnt - 1 && patchIdxJ == patchCnt - 1)
        diffPatchFilter = true;

    return diffPatchFilter;
}

bool mesh3D::filterPointTrig(const unsigned int surfVertIdx, const unsigned int surfTrigIdx) const
{
    int vI = this->surfVerts[surfVertIdx];
    const Vector4i& triI = this->surface[surfTrigIdx];

    bool adjacent = (vI == triI[0] || vI == triI[1] || vI == triI[2]);

    bool filtered = this->pointTrigFilterSet.contains(compressPair(surfVertIdx, surfTrigIdx));

    if (!this->characterAnimation)
        return (adjacent);

    unsigned int vertPatchIdx = 0, trigPatchIdx = 0;
    unsigned int patchCnt = this->patchSurfaceSep.size();
    for (unsigned int patchIdx = 0; patchIdx != patchCnt; ++patchIdx)
    {
        if (surfVertIdx >= this->patchVertexSep[patchIdx])
            vertPatchIdx = patchIdx;
        if (surfTrigIdx >= this->patchSurfaceSep[patchIdx])
            trigPatchIdx = patchIdx;
    }

    // specific for character garment animation!
    bool diffPatchFilter = AnimFilter(vertPatchIdx, trigPatchIdx, patchCnt);

    return (adjacent || diffPatchFilter);
}

using EdgePair = pair<uint64_t, uint64_t>;

bool mesh3D::filterEdgeEdge(const unsigned int surfEdgeIdxI, const unsigned int surfEdgeIdxJ) const
{
    const EdgePair& edgeI = this->surfEdges[surfEdgeIdxI];
    const EdgePair& edgeJ = this->surfEdges[surfEdgeIdxJ];

    bool adjacent = (edgeI.first == edgeJ.first || edgeI.first == edgeJ.second ||
                     edgeI.second == edgeJ.first || edgeI.second == edgeJ.second);
    bool filtered = 
        this->edgeEdgeFilterSet.contains(compressPair(surfEdgeIdxI, surfEdgeIdxJ)) ||
        this->edgeEdgeFilterSet.contains(compressPair(surfEdgeIdxJ, surfEdgeIdxI));

    if (!this->characterAnimation)
        return (adjacent);

    unsigned int edgeIPatchIdx = 0, edgeJPatchIdx = 0;
    unsigned int patchCnt = this->patchSurfaceSep.size();
    for (unsigned int patchIdx = 0; patchIdx != patchCnt; ++patchIdx)
    {
        if (surfEdgeIdxI >= this->patchEdgeSep[patchIdx])
            edgeIPatchIdx = patchIdx;
        if (surfEdgeIdxJ >= this->patchEdgeSep[patchIdx])
            edgeJPatchIdx = patchIdx;
    }

    // specific for character garment animation!
    bool diffPatchFilter = AnimFilter(edgeIPatchIdx, edgeJPatchIdx, patchCnt);

    return (adjacent || diffPatchFilter);
}

bool mesh3D::filterEdgeTrig(const unsigned int surfEdgeIdx, const unsigned int surfTrigIdx) const
{
    const EdgePair& edge = this->surfEdges[surfEdgeIdx];
    const Vector4i& triI = this->surface[surfTrigIdx];

    bool adjacent = (edge.first == triI[0] || edge.first == triI[1] || edge.first == triI[2] ||
        edge.second == triI[0] || edge.second == triI[1] || edge.second == triI[2]);
    bool filtered = this->edgeTrigFilterSet.contains(compressPair(surfEdgeIdx, surfTrigIdx));

    if (!this->characterAnimation)
        return (adjacent);

    unsigned int edgePatchIdx = 0, trigPatchIdx = 0;
    unsigned int patchCnt = this->patchSurfaceSep.size();
    for (unsigned int patchIdx = 0; patchIdx != patchCnt; ++patchIdx)
    {
        if (surfEdgeIdx >= this->patchEdgeSep[patchIdx])
            edgePatchIdx = patchIdx;
        if (surfTrigIdx >= this->patchSurfaceSep[patchIdx])
            trigPatchIdx = patchIdx;
    }

    // specific for character garment animation!
    bool diffPatchFilter = AnimFilter(edgePatchIdx, trigPatchIdx, patchCnt);

    return (adjacent || diffPatchFilter);
}

}