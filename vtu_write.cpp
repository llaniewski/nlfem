#include <stdio.h>
#include <stdlib.h>
#include <vector>

typedef std::vector<double> vec;


void write_vtu(char* filename, vec& points, vec& v, std::vector<size_t>& triangles) {
    FILE* f = fopen(filename, "w");
    if (f == NULL) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(-1);
    }
    size_t npoints = points.size()/3;
    size_t ncells = triangles.size()/6;
    if ((points.size() % 3 != 0) || (triangles.size() % 6 != 0) || (v.size() != points.size())) {
        fprintf(stderr, "Wrong number of points or triangles\n");
        exit(-1);
    }
    fprintf(f, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(f, "  <UnstructuredGrid>\n");
    fprintf(f, "    <Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n", npoints, ncells);
    fprintf(f, "      <PointData>\n");
    fprintf(f, "        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (size_t i=0; i<v.size(); i++) {
        fprintf(f, " %.15lg", v[i]);
        if (i % 6 == 0) fprintf(f, "\n");
    }
    fprintf(f, "\n");
    fprintf(f, "        </DataArray>\n");
    fprintf(f, "      </PointData>\n");
    fprintf(f, "      <CellData>\n");
    fprintf(f, "      </CellData>\n");
    fprintf(f, "      <Points>\n");
    fprintf(f, "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (size_t i=0; i<points.size(); i++) {
        fprintf(f, " %.15lg", points[i]);
        if (i % 6 == 0) fprintf(f, "\n");
    }
    fprintf(f, "\n");
    fprintf(f, "        </DataArray>\n");
    fprintf(f, "      </Points>\n");
    fprintf(f, "      <Cells>\n");
    fprintf(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for (size_t i=0; i<triangles.size(); i++) {
        fprintf(f, " %ld", triangles[i]);
        if (i % 6 == 0) fprintf(f, "\n");
    }
    fprintf(f, "\n");
    fprintf(f, "        </DataArray>\n");
    fprintf(f, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for (size_t i=0; i<triangles.size(); i++) {
        fprintf(f, " %ld", (i+1)*6);
        if (i % 6 == 0) fprintf(f, "\n");
    }
    fprintf(f, "\n");
    fprintf(f, "        </DataArray>\n");
    fprintf(f, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (size_t i=0; i<triangles.size(); i++) {
        fprintf(f, " %d", 13);
        if (i % 6 == 0) fprintf(f, "\n");
    }
    fprintf(f, "\n");
    fprintf(f, "        </DataArray>\n");
    fprintf(f, "      </Cells>\n");
    fprintf(f, "    </Piece>\n");
    fprintf(f, "  </UnstructuredGrid>\n");
    fprintf(f, "</VTKFile>\n");
    fprintf(f, "\n");
    fclose(f);
}