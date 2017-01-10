//editor.execute( new AddObjectCommand( mesh ) );
var nDefVols = 0;
var matList = [];
var nDefMats = 0;
var colors = [0xFFFFFF, 0x000000, 0x333333, 0x666666, 0x999999, 0xCCCCCC, 0xCCCC99, 0x9999CC, 0x666699, 0x660000, 0x663300, 0x996633, 0x003300, 0x003333, 0x003399, 0x000066, 0x330066, 0x660066, 0x990000, 0x993300, 0xCC9900, 0x006600, 0x336666, 0x0033FF, 0x000099, 0x660099, 0x990066, 0xCC0000, 0xCC3300, 0xFFCC00, 0x009900, 0x006666, 0x0066FF, 0x0000CC, 0x663399, 0xCC0099, 0xFF0000, 0xFF3300, 0xFFFF00, 0x00CC00, 0x009999, 0x0099FF, 0x0000FF, 0x9900CC, 0xFF0099, 0xCC3333, 0xFF6600, 0xFFFF33, 0x00FF00, 0x00CCCC, 0x00CCFF, 0x3366FF, 0x9933FF, 0xFF00FF, 0xFF6666, 0xFF6633, 0xFFFF66, 0x66FF66, 0x66CCCC, 0x00FFFF, 0x3399FF, 0x9966FF, 0xFF66FF, 0xFF9999, 0xFF9966, 0xFFFF99, 0x99FF99, 0x66FFCC, 0x99FFFF, 0x66CCFF, 0x9999FF, 0xFF99FF, 0xFFCCCC, 0xFFCC99, 0xFFFFCC, 0xCCFFCC, 0x99FFCC, 0xCCFFFF, 0x99CCFF, 0xCCCCFF, 0xFFCCFF];
var colorNames = ["color1", "color2", "color3", "color4", "color5", "color6", "color7", "color8", "color9", "color10", "color11", "color12", "color13", "color14", "color15", "color16", "color17", "color18", "color19", "color20", "color21", "color22", "color23", "color24", "color25", "color26", "color27", "color28", "color29", "color30", "color31", "color32", "color33", "color34", "color35", "color36", "color37", "color38", "color39", "color40", "color41", "color42", "color43", "color44", "color45", "color46", "color47", "color48", "color49", "color50", "color51", "color52", "color53", "color54", "color55", "color56", "color57", "color58", "color59", "color60", "color61", "color62", "color63", "color64", "color65", "color66", "color67", "color68", "color69", "color70", "color71", "color72", "color73", "color74", "color75", "color76", "color77", "color78", "color79", "color80", "color81"];
var nColors = 81;

var cloneIndex = [];
var cloneCopy = [];
var cloneMother = [];
var cloneCopyName = [];
var cloneV = [[]];
var nClonedVolumes = 0;
var cloneMerge = [];

//Mesh precision values
var narcseg=20;
var ntubseg=10;
var ncseg=45; //number of segments on the circonference
var nzseg=45; //number of segments on Z
var nphiseg=45;
var ntheseg=45;


// Geolist will contail all geometries.
var geoList = [];
// Geolist will contail all geometries.
var meshList = [];

// Alll the tree structure of physical volumes will has "root" as the root of the tree
// Any operation will be performed using three-js functions.
var root = new THREE.Object3D();
var root.name = "root";

var spinner;

var clipPlanes = [
    new THREE.Plane(new THREE.Vector3(-1, 0, 0), 0),
    new THREE.Plane(new THREE.Vector3(0, -1, 0), 0),
    new THREE.Plane(new THREE.Vector3(0, 0, -1), 0)
];

function setColors(c) {
    for (i = 0; i < c.length; i++) {
        colors[i] = c[i].cl;
        colorNames[i] = c[i].cn;
    }
    nColors = colors.length;
    return nColors;
}

function addNewColor(c, n) {
    colorIndex = -1;
    for (i = 0; i < nColors; i++) {
        if (c == colors[i]) {
            colorIndex = i;
            break;
        }
    }
    if (colorIndex < 0) {
        colorIndex = nColors;
        colorNames[nColors] = n;
        colors[nColors++] = c;
       defColorMat(colorIndex);
    }
    return colorIndex;
}


function getMat(index) {
    if (index < nDefMats) {
        return matList[index];
    } else {
        //console.log("***Warning Material undefined material index:"+index+" chosen random material");
        index = Math.floor(nDefMats * Math.random());
        return matList[index];
    }

}

function defColorMat(i=0,transparency = false, opacity = 0.6, side = THREE.DoubleSide) {
        matList[nDefMats++] = new THREE.MeshPhongMaterial({
            name: colorNames[i],
            color: colors[i],
            transparent: transparency,
            opacity: opacity,
            side: side,
            shininess: 100,
            clippingPlanes: clipPlanes,
            clipIntersection: true,
            clipShadows: true
        });

}


function defColorMaterials(transparency = false, opacity = 0.6, side = THREE.DoubleSide) {
    for (i = 0; i < nColors; i++) {
        defColorMat(i,transparency, opacity, side) {
    }
    return nDefMats;
}



function addMeshToList(mesh) {
    meshList[nDefVols++] = mesh;
    return nDefVols;
}


function makeClone(motherName, volumeName, copyNumber, v) {
    var vId = root.getObjectByName(volumeName);
    if (!vId) {
        console.log("*****ERROR Clone volume:" + name + " not defined ...");
        return 0;
    }
    var mesh = vId.clone();
    setMeshMatrix(mesh, v);

    mesh.name = volumeName + "-" + copyNumber;

    getMother(motherName).add(mesh);

    addMeshToList(mesh);

    return mesh;

}


function getMother(volName) {
    if (volName == "root") {
        mother = root;
    } else {
        mother = root.getObjectByName(volName);
    }
    if (mother) {
        return mother;
    } else {
        return root;
    }
}

function setMeshMatrix(mesh, v) {
    if (v.length > 0) {
        var position = new THREE.Vector3();
        var quaternion = new THREE.Quaternion();
        var scale = new THREE.Vector3();
        var m = new THREE.Matrix4().set(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[14], v[14], v[15]);
        m.decompose(position, quaternion, scale);
        mesh.position.copy(position);
        mesh.scale.copy(scale);
        mesh.quaternion.copy(quaternion);
    }
}

function setMeshMatrixC(mesh, v) {
    if (v.length > 0) {
        var position = new THREE.Vector3();
        var quaternion = new THREE.Quaternion();
        var scale = new THREE.Vector3();
        var m = new THREE.Matrix4().set(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[14], v[14], v[15]);
        m.decompose(position, quaternion, scale);
        mesh.position.copy(position);
        mesh.scale.copy(scale);
        //mesh.quaternion.copy(quaternion);
    }
}

function setMeshProperties(mesh, mother,name,v,visible) {
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    mesh.visible = visible;
    mesh.name = name;

    setMeshMatrix(mesh, v);

    getMother(mother).add(mesh);
    addMeshToList(mesh);
}

function makeComposite(mother, mat, visible, v, Vol1, Vol2, operation, newVol) {
    compositeL = root.getObjectByName(Vol1);
    compositeR = root.getObjectByName(Vol2);

    if (operation == "Union1") {

        var geometry = new THREE.Geometry();
        compositeL.geometry;
        geometry.merge(compositeL.geometry, compositeL.matrix);
        geometry.merge(compositeR.geometry, compositeR.matrix);
        var mesh = new THREE.Mesh(geometry, getMat(mat));

    } else {
        CleftBSP = new ThreeBSP(compositeL);
        CrightBSP = new ThreeBSP(compositeR);
        if (operation == "Intersection") {
            resultBSP = CleftBSP.intersect(CrightBSP);

        } else if (operation == "Union") {
            resultBSP = CleftBSP.union(CrightBSP);

        } else if (operation == "Subtraction") {
            resultBSP = CleftBSP.subtract(CrightBSP);
        }
        var mesh = resultBSP.toMesh();
    }

    mesh.name = newVol
    mesh.material = compositeL.material;
    mesh.geometry.computeFaceNormals();
    mesh.geometry.computeVertexNormals();


    mesh.castShadow = true;
    mesh.receiveShadow = true;
    mesh.visible = visible;
    setMeshMatrixC(mesh, v);
    var motherComp = getMother(mother);
    motherComp.add(mesh);
    addMeshToList(mesh);

    motherComp.remove(compositeL);
    motherComp.remove(compositeR);

    return mesh;
}

function geoCreateMesh(geoId, mother, mat, visible, name, v) {
    var mesh = new THREE.Mesh(geoList[geoId], getMat(mat));
    setMeshProperties(mesh, mother,name,v,visible);
    return mesh;
}




function geoCreateBox(geoId, x, y, z) {
    geoList[geoId] = new THREE.BoxGeometry(x, y, z);
    return geoList[geoId];
}



function geoCreateLathe(geoId, r, z, phiStart = 0, phiLength = 2 * Math.PI) {
    if ((Math.PI * 2 - phiLength) < 0.0001) phiLength = Math.PI * 2;
    var points = [];
    var np = r.length;
    for (i = 0; i < np; i++) {
        points.push(new THREE.Vector2(r[i], z[i]));
        //console.log(r[i],z[i]);
    }

    // Create the geometry & mesh in PrimFile
    var geom = new THREE.LatheGeometry(points, ncseg, phiStart, phiLength);

    if (phiLength < 2 * Math.PI) {
        var myShape = new THREE.Shape();
        myShape.moveTo(r[0], z[0]);

        for (var i = 1; i < np; i++) myShape.lineTo(r[i], z[i]);
        var sgeom = new THREE.ShapeGeometry(myShape);
        sgeom.rotateY(-0.5 * Math.PI + phiStart);
        geom.merge(sgeom);
        sgeom.rotateY(phiLength);
        geom.merge(sgeom);
    }

    // rotate to position it correctly along Z axis (in THREE-JS the parabloid is created along the Y axis)
    geom.rotateX(0.5 * Math.PI);
    geom.rotateZ(0.5 * Math.PI);

    geoList[geoId] = geom;
    return geoList[geoId];
}


function geoCreateTorus(geoId, R, Rmin, Rmax, phiStart = 0, phiLength = 2 * Math.pi) {

    if ((Math.PI * 2 - phiLength) < 0.0001) phiLength = Math.PI * 2;
    var geom = new THREE.TorusGeometry(R, Rmax, ntubseg, narcseg, phiLength);

    if (phiLength < 2 * Math.pi) {

        var m = new THREE.Matrix4();
        var geom1 = new THREE.TorusGeometry(R, Rmin, ntubseg, narcseg, phiLength);
        geom.merge(geom1);

        m.makeRotationX(phiLength);
        m.setPosition(new THREE.Vector3(R, 0, 0));

        var ring = new THREE.RingGeometry(Rmin, Rmax, ntubseg, 8);
        geom.merge(ring, m);

        m.identity();

        m.makeRotationX(Math.pi / 2);
        m.makeRotationY(phiLength);
        m.setPosition(new THREE.Vector3(R * Math.cos(phiLength), R * Math.sin(phiLength), 0));
        geom.merge(ring, m);
    }

    geom.rotateZ(phiStart);

    if (phiLength < 2 * Math.pi) {

        var m = new THREE.Matrix4();
        var geom1 = new THREE.TorusGeometry(R, Rmin, ntubseg, narcseg, phiLength);
        geom.merge(geom1);

        m.makeRotationX(phiLength);
        m.setPosition(new THREE.Vector3(R, 0, 0));

        var ring = new THREE.RingGeometry(Rmin, Rmax, ntubseg, 8);
        geom.merge(ring, m);

        m.identity();

        m.makeRotationX(Math.pi / 2);
        m.makeRotationY(phiLength);
        m.setPosition(new THREE.Vector3(R * Math.cos(phiLength), R * Math.sin(phiLength), 0));
        geom.merge(ring, m);
    }
    geom.rotateZ(phiStart);

    geoList[geoId] = geom;
    return geoList[geoId];

}

function geoCreateExtruded(geoId,name, x, y, z, px, py, f) {
    np = x.length - 1;
    nz = z.length;
    //define the top and bottom faces to be added to the geometry
    var shape = new THREE.Shape();
    var shape1 = new THREE.Shape();
    shape.moveTo(x[0] * f[0] + px[0], y[0] * f[0] + py[0])




    shape1.moveTo(x[0] * f[nz - 1] + px[nz - 1], y[0] * f[nz - 1] + py[nz - 1])
    for (i = 1; i <= np; i++) {
        shape.lineTo(x[i] * f[0] + px[0], y[i] * f[0] + py[0]);
        shape1.lineTo(x[i] * f[nz - 1] + px[nz - 1], y[i] * f[nz - 1] + py[nz - 1])
    }
    if ((nz == 2) && (f[0] * f[1] == 1)) {
        var extrudeSettings = {
            amount: z[1] - z[0],
            steps: 1,
        };

        var geom = new THREE.ExtrudeGeometry(shape, extrudeSettings);

    } else {

        var shapeGeom = new THREE.ShapeGeometry(shape);
        var shapeGeom1 = new THREE.ShapeGeometry(shape1);
        shapeGeom.translate(0, 0, z[0]); // bottom shape
        shapeGeom1.translate(0, 0, z[nz - 1]); // top shape

        // define the lateral shapes defined by the blue-print coordinates 
        var geom = new THREE.Geometry();
        for (j = 0; j < nz; j++) {
            for (i = 0; i < np; i++) {
                geom.vertices.push(new THREE.Vector3(x[i] * f[j] + px[j], y[i] * f[j] + py[j], z[j]));
            }
        }
        for (j = 0; j < nz - 1; j++) {
            for (i = 0; i < np - 1; i++) {
                geom.faces.push(new THREE.Face3(j * np + i, j * np + i + 1, (j + 1) * np + i + 1));
                geom.faces.push(new THREE.Face3((j + 1) * np + i + 1, (j + 1) * np + i, j * np + i));
            }
        }
        for (j = 1; j < nz; j++) {
            geom.faces.push(new THREE.Face3(j * np - 1, j * np - np, j * np));
            geom.faces.push(new THREE.Face3(j * np, (j + 1) * np - 1, j * np - 1));
        }
        // merge top and bottom shapes with the lateral shape to compose the ful volume

        geom.merge(shapeGeom);
        geom.merge(shapeGeom1);
        geom.computeFaceNormals();

    }

    geoList[geoId] = geom;
    return geoList[geoId];
}


function geoCreateSphere(geoId, name, Rmin, Rmax, thetaStart = 0, thetaLength = 2 * Math.PI, phiStart = 0, phiLength = 2 * Math.PI) {
    if ((Math.PI * 2 - phiLength) < 0.0001) phiLength = Math.PI * 2;
    if ((Math.PI - thetaLength) < 0.0001) thetaLength = Math.PI;
    if ((thetaLength == Math.PI) && (phiLength == 2 * Math.PI)) {
        var geom = new THREE.SphereGeometry(Rmax, nphiseg, ntheseg);
        if (Rmin < Rmax) {
            var mesh1 = new THREE.Mesh(geom);
            var mesh2 = new THREE.Mesh(new THREE.SphereGeometry(Rmin, nphiseg, ntheseg));
            var CleftBSP = new ThreeBSP(mesh1);
            var CrightBSP = new ThreeBSP(mesh2);
            var resultBSP = CleftBSP.subtract(CrightBSP);
            var mesh = resultBSP.toMesh();
            geom = mesh.geometry;
        }

    } else {

        var r = [];
        var z = [];

        r[0] = Rmax * Math.sin(thetaStart);
        z[0] = Rmax * Math.cos(thetaStart);
        r[ntheseg] = Rmax * Math.sin(thetaStart + thetaLength);
        z[ntheseg] = Rmax * Math.cos(thetaStart + thetaLength);

        var td, theta;
        td = thetaLength / ntheseg;

        theta = thetaStart;

        for (nseg = 1; nseg < ntheseg; nseg++) {
            theta += td;
            r[nseg] = Rmax * Math.sin(theta);
            z[nseg] = Rmax * Math.cos(theta);
        }

        if (Rmin > 0) {

            theta1 = thetaStart + thetaLength;
            nt = ntheseg + 1;
            for (nseg = 1; nseg < ntheseg; nseg++) {
                nt++;
                theta1 -= td;
                r[nt] = Rmin * Math.sin(theta1);
                z[nt] = Rmin * Math.cos(theta1);
            }
            r[ntheseg + 1] = Rmin * Math.sin(thetaStart + thetaLength);
            z[ntheseg + 1] = Rmin * Math.cos(thetaStart + thetaLength);

            r[nt + 1] = Rmin * Math.sin(thetaStart);
            z[nt + 1] = Rmin * Math.cos(thetaStart);

            r[nt + 2] = r[0];
            z[nt + 2] = z[0];

        } else {

            r[ntheseg + 1] = 0;
            z[ntheseg + 1] = Rmax * Math.cos(thetaStart + thetaLength);

            r[ntheseg + 2] = r[0];
            z[ntheseg + 2] = z[0];
        }

        var points = [];
        np = r.length;
        for (i = 0; i < np; i++) points.push(new THREE.Vector2(r[i], z[i]));

        // Create the geometry & mesh in PrimFile
        var geom = new THREE.LatheGeometry(points, nphiseg, phiStart, phiLength);

        if (phiLength < 2 * Math.PI) {
            var myShape = new THREE.Shape();
            myShape.moveTo(r[0], z[0]);
            //console.log(r[0],z[0]);
            for (i = 1; i < np; i++) {
                myShape.lineTo(r[i], z[i]);
                //console.log(r[i],z[i]);
            }
            var sgeom = new THREE.ShapeGeometry(myShape);
            sgeom.rotateY(-0.5 * Math.PI + phiStart);
            geom.merge(sgeom);
            sgeom.rotateY(phiLength);
            geom.merge(sgeom);
        }
    }
    geom.rotateX(Math.PI / 2);
    geom.rotateZ(Math.PI / 2);
    geom.computeFaceNormals();

    geoList[geoId] = geom;
    return geoList[geoId];
}


function geoCreateParabloid(geoId, Rlo, Rhi, Dz) {
    var r = [];
    var z = [];
    //start from the center to include lower & upper surface at -Dz & Dz
    r[0] = 0;
    z[0] = -Dz;
    var countLow = 0;
    if (Rlo > 0.) {
        countLow = 1;
        r[1] = Rlo;
        z[1] = -Dz;
    }

    r[nzSegments + 1] = Rhi;
    z[nzSegments + 1] = Dz;
    r[nzSegments + 2] = 0;
    z[nzSegments + 2] = Dz;

    // Calculate A and B constants from parabloid equations: (-Dz = A*Rlo*Rlo + B) & (Dz = A*Rhi*Rhi + B)
    var d2 = (Rhi * Rhi - Rlo * Rlo);
    var dd = 1.0 / d2;
    var A = 2.0 * Dz * dd;
    var B = -Dz * (Rlo * Rlo + Rhi * Rhi) * dd;

    // For a given value of Z between [-Dz,Dz], calculate R  
    // from parabloid equations: (-Dz = a*Rlo*Rlo + b) & (Dz = a*Rhi*Rhi + b)
    // the surface of revolution of a parabola described by: `z = a*(x*x + y*y) + b`
    // We need only the points of x or y axis which is the radious ar a pont z and we name it (R)
    // R= SQRT((Z-B)/A).

    var RR, ZZ;

    for (nseg = 1; nseg < nzSegments; nseg++) {
        ZZ = -Dz + nseg * 2 * Dz / nzSegments;
        RR = Math.sqrt((ZZ - B) / A);
        r[nseg + countLow] = RR;
        z[nseg + countLow] = ZZ;
    }
    return geoCreateLathe(geoId, r, z );

}


function geoCreateHyperbloid(geoId, Rmin, Rmax, StIn, StOut, Dz) {
    var r = [];
    var z = [];
    // Calculate Rin and Rout for each Z position from Hyperbloid equations: 
    // Rin=SQRT( Rmin**2 + (tan(Stin)*Z)**2 )
    // Rout=SQRT( Rmax**2 + (tan(Stout)*Z)**2 )

    var ZZ = -Dz;
    var sd = 2 * Dz / nzSegments;
    linebk = nzSegments + 1;
    for (nseg = 0; nseg <= nzSegments; nseg++) {
        var PIZ = Math.tan(StIn) * ZZ;
        var POZ = Math.tan(StOut) * ZZ;
        var Rin = Math.sqrt(Math.pow(Rmin, 2.) + Math.pow(PIZ, 2.));
        var Rout = Math.sqrt(Math.pow(Rmax, 2.) + Math.pow(POZ, 2.));
        r[nseg] = Rout;
        z[nseg] = ZZ;
        r[2 * nzSegments + 1 - nseg] = Rin;
        z[2 * nzSegments + 1 - nseg] = ZZ;

        ZZ += sd;
    }

    r[2 * nzSegments + 2] = r[0];
    z[2 * nzSegments + 2] = z[0];

    return return geoCreateLathe(geoId, r, z );
}


function geoCreateTrapezoid(geoId, x, y, z, nz) {
    // Define a generic paralelepiped, trapezoid etc
    //  x,y,z: coordinates of bottom [index:0-3] and top  [index:4-7] faces
    // nz=number of divisions in Z
    // The resulting volume is a set of vertexes and faces
    //

    if (nz == 1) nz = 2;
    var geom = new THREE.Geometry();
    var cur = [];
    var tcp = [];
    for (i = 0; i < 4; i++) {
        cur[i] = new THREE.LineCurve3(new THREE.Vector3(x[i], y[i], z[i]), new THREE.Vector3(x[i + 4], y[i + 4], z[i + 4]));
        tcp[i] = cur[i].getPoints(nz - 1);
    }

    for (j = 0; j < nz; j++) {
        for (i = 0; i < 4; i++) {
            geom.vertices.push(tcp[i][j]);
        }
    }
    geom.faces.push(new THREE.Face3(0, 1, 2));
    geom.faces.push(new THREE.Face3(2, 3, 0));
    geom.faces.push(new THREE.Face3(4 * (nz - 1), 4 * (nz - 1) + 1, 4 * (nz - 1) + 2));
    geom.faces.push(new THREE.Face3(4 * (nz - 1) + 2, 4 * (nz - 1) + 3, 4 * (nz - 1)));
    var f4 = [0, 0, 0, 4];
    np = 4 * (nz - 1);
    for (i = 0; i < 4; i++) {
        for (j = 0; j < np; j += 4) {
            A = i + j;
            B = A + 1;
            C = B + 4;
            geom.faces.push(new THREE.Face3(A, B - f4[i], C - f4[i]));
            geom.faces.push(new THREE.Face3(C - f4[i], C - 1, A));
        }
    }
    geom.computeFaceNormals();
    geom.computeVertexNormals();

    geoList[geoId] = geom;
    return geoList[geoId];
}


function geoCreatePara(geoId, dx, dy, dz, alpha, theta, phi) {
    var fxy = Math.tan(alpha);
    var tth = Math.tan(theta);
    var fxz = tth * Math.cos(phi);
    var fyz = tth * Math.sin(phi);
    var dx1 = dx + dy * Math.abs(fxy) + dz * Math.abs(fxz);
    var dy1 = dy + dz * Math.abs(fyz);
    var dz1 = dz;
    var x = [-dx1, -dx1 + 2 * dx, -dx1 + 2 * dx + fxy * 2 * dy, -dx1 + fxy * 2 * dy, dx1 - 2 * dx - fxy * 2 * dy, dx1 - fxy * 2 * dy, dx1, dx1 - fxy * 2 * dy];
    var y = [-dy1, -dy1, -dy1 + 2 * dy, -dy1 + 2 * dy, dy1 - 2 * dy, dy1 - 2 * dy, dy1, dy1];
    var z = [-dz, -dz, -dz, -dz, dz, dz, dz, dz];

    return geoCreateTrapezoid(geoId, x, y, z, 8);

}

function geoCreateTrap(geoId, dz, theta, phi, bl1, tl1, bl2, tl2, h1, h2, alpha1, alpha2) {
    var tx = Math.tan(theta) * Math.cos(phi);
    var ty = Math.tan(theta) * Math.sin(phi);
    var ta1 = Math.tan(alpha1);
    var ta2 = Math.tan(alpha2);

    var x = [-dz * tx - h1 * ta1 - bl1, -dz * tx + h1 * ta1 - tl1, -dz * tx + h1 * ta1 + tl1, -dz * tx - h1 * ta1 + bl1, dz * tx - h2 * ta2 - bl2, dz * tx + h2 * ta2 - tl2, dz * tx + h2 * ta2 + tl2, dz * tx - h2 * ta2 + bl2];
    var y = [-dz * ty - h1, -dz * ty + h1, -dz * ty + h1, -dz * ty - h1, dz * ty - h2, dz * ty + h2, dz * ty + h2, dz * ty - h2];
    var z = [-dz, -dz, -dz, -dz, dz, dz, dz, dz];

    return geoCreateTrapezoid(geoId, x, y, z, 8);;
}



function geoCreateGTra(geoId, dz, theta, phi, bl1, tl1, bl2, tl2, h1, h2, alpha1, alpha2, twist) {

    //same like TGeoTrap - begin
    var tx = Math.tan(theta) * Math.cos(phi);
    var ty = Math.tan(theta) * Math.sin(phi);
    var ta1 = Math.tan(alpha1);
    var ta2 = Math.tan(alpha2);

    var x = [-dz * tx - h1 * ta1 - bl1, -dz * tx + h1 * ta1 - tl1, -dz * tx + h1 * ta1 + tl1, -dz * tx - h1 * ta1 + bl1, dz * tx - h2 * ta2 - bl2, dz * tx + h2 * ta2 - tl2, dz * tx + h2 * ta2 + tl2, dz * tx - h2 * ta2 + bl2];
    var y = [-dz * ty - h1, -dz * ty + h1, -dz * ty + h1, -dz * ty - h1, dz * ty - h2, dz * ty + h2, dz * ty + h2, dz * ty - h2];
    var z = [-dz, -dz, -dz, -dz, dz, dz, dz, dz];

    //same like TGeoTrap - end

    var x0, y0;
    var th = theta;
    var ph = phi;
    // Coordinates of the center of the bottom face
    var xc = -dz * Math.sin(th) * Math.cos(ph);
    var yc = -dz * Math.sin(th) * Math.sin(ph);

    for (i = 0; i < 4; i++) {
        x0 = x[i] - xc;
        y0 = y[i] - yc;
        x[i] = x0 * Math.cos(-0.5 * twist) + y0 * Math.sin(-0.5 * twist) + xc;
        y[i] = -x0 * Math.sin(-0.5 * twist) + y0 * Math.cos(-0.5 * twist) + yc;
    }
    // Coordinates of the center of the top face
    xc = -xc;
    yc = -yc;
    for (i = 4; i < 8; i++) {
        x0 = x[i] - xc;
        y0 = y[i] - yc;
        x[i] = x0 * Math.cos(0.5 * twist) + y0 * Math.sin(0.5 * twist) + xc;
        y[i] = -x0 * Math.sin(0.5 * twist) + y0 * Math.cos(0.5 * twist) + yc;
    }
    return geoCreateTrapezoid(geoId, x, y, z, 8);;
}


function geoCreateTrd(geoId, dx1, dx2, dy1, dy2, dz) {
    var x = [-dx1, dx1, dx1, -dx1, -dx2, dx2, dx2, -dx2];
    var y = [dy1, dy1, -dy1, -dy1, dy2, dy2, -dy2, -dy2];
    var z = [-dz, -dz, -dz, -dz, dz, dz, dz, dz];
    return     return geoCreateTrapezoid(geoId, x, y, z, 8);;
}


function geoCreateConeSeg(geoId, Dz, Rmin1, Rmin2, Rmax1, Rmax2, phiStart = 0, phiLength = 2 * Math.PI) {
    if ((Math.PI * 2 - phiLength) < 0.0001) {
        phiLength = Math.PI * 2;
    }
    //console.log(Rmin1,Rmax1,Rmax2,Rmin2,Rmin1);
    var r = [Rmin1, Rmax1, Rmax2, Rmin2, Rmin1];
    var z = [-Dz, -Dz, Dz, Dz, -Dz];

    return geoCreateLathe(geoId, r, z, phiStart, phiLength);
}


function geoCreateCylinder(geoId, Dz, Rtop, Rbottom, phiStart = 0, phiLength = 2 * Math.PI) {
    if ((Math.PI * 2 - phiLength) < 0.0001) phiLength = Math.PI * 2;

    var geom = new THREE.CylinderGeometry(Rtop, Rbottom, Dz, ncseg, nzseg, true, phiStart, phiLength);
    geom.rotateX(0.5 * Math.PI);
    var kp = geoList.length;
    geoList[geoId] = geom;
    return geoList[geoId];
}


function geoCreateElipticalCylinder(geoId, Dz, A, B, phiStart = 0, phiLength = 2 * Math.PI) {
    if ((Math.PI * 2 - phiLength) < 0.0001) phiLength = Math.PI * 2;
    if (geom == null) geom = new THREE.CylinderGeometry(1, 1, Dz, ncSegments, nzSegments, true, phiStart, phiLength);
    geom.scale(A, 1, B);
    geom.rotateX(0.5 * Math.PI);
    geoList[geoId] = geom;
    return geoList[geoId];
}

function drawElipticalCylinder(geoId, mother, mat, visible, name, v) {
    if (geom == null) geom = geoCreateElipticalCylinder(geom, name, v, ncSegments, nzSegments, Dz, A, B, phiStart, phiLength);
    var mesh = new THREE.Mesh(geom, getMat(mat));
    setMeshProperties(mesh, mother,name,v,visible);

    return mesh;
}

function drawGroup(geoId, mother, visible, name, v) {
    var group = new THREE.Object3D()
    group.name = name;
    group.visible = visible;
    setMeshMatrix(group, v);
    getMother(mother).add(group);
    addMeshToList(group);
    return group;
}

function showSpinner() {

    var opts = {
        lines: 13, // The number of lines to draw
        length: 20, // The length of each line
        width: 10, // The line thickness
        radius: 30, // The radius of the inner circle
        corners: 1, // Corner roundness (0..1)
        rotate: 0, // The rotation offset
        direction: 1, // 1: clockwise, -1: counterclockwise
        color: '#ffff00', // #rgb or #rrggbb or array of colors
        speed: 1, // Rounds per second
        trail: 60, // Afterglow percentage
        shadow: false, // Whether to render a shadow
        hwaccel: false, // Whether to use hardware acceleration
        className: 'spinner', // The CSS class to assign to the spinner
        zIndex: 2e9, // The z-index (defaults to 2000000000)
        top: 'auto', // Top position relative to parent in px
        left: 'auto' // Left position relative to parent in px
    };
    var target = document.getElementById('WebGL-output');
    spinner = new Spinner(opts).spin(target);
    return spinner;
}

function hideSpinner(spinner) {
    spinner.stop();
}


(function(t, e) {
    if (typeof exports == "object") module.exports = e();
    else if (typeof define == "function" && define.amd) define(e);
    else t.Spinner = e()
})(this, function() {
    "use strict";
    var t = ["webkit", "Moz", "ms", "O"],
        e = {},
        i;

    function o(t, e) {
        var i = document.createElement(t || "div"),
            o;
        for (o in e) i[o] = e[o];
        return i
    }

    function n(t) {
        for (var e = 1, i = arguments.length; e < i; e++) t.appendChild(arguments[e]);
        return t
    }
    var r = function() {
        var t = o("style", { type: "text/css" });
        n(document.getElementsByTagName("head")[0], t);
        return t.sheet || t.styleSheet
    }();

    function s(t, o, n, s) {
        var a = ["opacity", o, ~~(t * 100), n, s].join("-"),
            f = .01 + n / s * 100,
            l = Math.max(1 - (1 - t) / o * (100 - f), t),
            u = i.substring(0, i.indexOf("Animation")).toLowerCase(),
            d = u && "-" + u + "-" || "";
        if (!e[a]) {
            r.insertRule("@" + d + "keyframes " + a + "{" + "0%{opacity:" + l + "}" + f + "%{opacity:" + t + "}" + (f + .01) + "%{opacity:1}" + (f + o) % 100 + "%{opacity:" + t + "}" + "100%{opacity:" + l + "}" + "}", r.cssRules.length);
            e[a] = 1
        }
        return a
    }

    function a(e, i) {
        var o = e.style,
            n, r;
        i = i.charAt(0).toUpperCase() + i.slice(1);
        for (r = 0; r < t.length; r++) {
            n = t[r] + i;
            if (o[n] !== undefined) return n
        }
        if (o[i] !== undefined) return i
    }

    function f(t, e) {
        for (var i in e) t.style[a(t, i) || i] = e[i];
        return t
    }

    function l(t) {
        for (var e = 1; e < arguments.length; e++) {
            var i = arguments[e];
            for (var o in i)
                if (t[o] === undefined) t[o] = i[o]
        }
        return t
    }

    function u(t) {
        var e = { x: t.offsetLeft, y: t.offsetTop };
        while (t = t.offsetParent) e.x += t.offsetLeft, e.y += t.offsetTop;
        return e
    }

    function d(t, e) {
        return typeof t == "string" ? t : t[e % t.length]
    }
    var p = { lines: 12, length: 7, width: 5, radius: 10, rotate: 0, corners: 1, color: "#000", direction: 1, speed: 1, trail: 100, opacity: 1 / 4, fps: 20, zIndex: 2e9, className: "spinner", top: "auto", left: "auto", position: "relative" };

    function c(t) {
        if (typeof this == "undefined") return new c(t);
        this.opts = l(t || {}, c.defaults, p)
    }
    c.defaults = {};
    l(c.prototype, {
        spin: function(t) {
            this.stop();
            var e = this,
                n = e.opts,
                r = e.el = f(o(0, { className: n.className }), { position: n.position, width: 0, zIndex: n.zIndex }),
                s = n.radius + n.length + n.width,
                a, l;
            if (t) {
                t.insertBefore(r, t.firstChild || null);
                l = u(t);
                a = u(r);
                f(r, { left: (n.left == "auto" ? l.x - a.x + (t.offsetWidth >> 1) : parseInt(n.left, 10) + s) + "px", top: (n.top == "auto" ? l.y - a.y + (t.offsetHeight >> 1) : parseInt(n.top, 10) + s) + "px" })
            }
            r.setAttribute("role", "progressbar");
            e.lines(r, e.opts);
            if (!i) {
                var d = 0,
                    p = (n.lines - 1) * (1 - n.direction) / 2,
                    c, h = n.fps,
                    m = h / n.speed,
                    y = (1 - n.opacity) / (m * n.trail / 100),
                    g = m / n.lines;
                (function v() {
                    d++;
                    for (var t = 0; t < n.lines; t++) {
                        c = Math.max(1 - (d + (n.lines - t) * g) % m * y, n.opacity);
                        e.opacity(r, t * n.direction + p, c, n)
                    }
                    e.timeout = e.el && setTimeout(v, ~~(1e3 / h))
                })()
            }
            return e
        },
        stop: function() {
            var t = this.el;
            if (t) {
                clearTimeout(this.timeout);
                if (t.parentNode) t.parentNode.removeChild(t);
                this.el = undefined
            }
            return this
        },
        lines: function(t, e) {
            var r = 0,
                a = (e.lines - 1) * (1 - e.direction) / 2,
                l;

            function u(t, i) {
                return f(o(), { position: "absolute", width: e.length + e.width + "px", height: e.width + "px", background: t, boxShadow: i, transformOrigin: "left", transform: "rotate(" + ~~(360 / e.lines * r + e.rotate) + "deg) translate(" + e.radius + "px" + ",0)", borderRadius: (e.corners * e.width >> 1) + "px" })
            }
            for (; r < e.lines; r++) {
                l = f(o(), { position: "absolute", top: 1 + ~(e.width / 2) + "px", transform: e.hwaccel ? "translate3d(0,0,0)" : "", opacity: e.opacity, animation: i && s(e.opacity, e.trail, a + r * e.direction, e.lines) + " " + 1 / e.speed + "s linear infinite" });
                if (e.shadow) n(l, f(u("#000", "0 0 4px " + "#000"), { top: 2 + "px" }));
                n(t, n(l, u(d(e.color, r), "0 0 1px rgba(0,0,0,.1)")))
            }
            return t
        },
        opacity: function(t, e, i) {
            if (e < t.childNodes.length) t.childNodes[e].style.opacity = i
        }
    });

    function h() {
        function t(t, e) {
            return o("<" + t + ' xmlns="urn:schemas-microsoft.com:vml" class="spin-vml">', e)
        }
        r.addRule(".spin-vml", "behavior:url(#default#VML)");
        c.prototype.lines = function(e, i) {
            var o = i.length + i.width,
                r = 2 * o;

            function s() {
                return f(t("group", { coordsize: r + " " + r, coordorigin: -o + " " + -o }), { width: r, height: r })
            }
            var a = -(i.width + i.length) * 2 + "px",
                l = f(s(), { position: "absolute", top: a, left: a }),
                u;

            function p(e, r, a) { n(l, n(f(s(), { rotation: 360 / i.lines * e + "deg", left: ~~r }), n(f(t("roundrect", { arcsize: i.corners }), { width: o, height: i.width, left: i.radius, top: -i.width >> 1, filter: a }), t("fill", { color: d(i.color, e), opacity: i.opacity }), t("stroke", { opacity: 0 })))) }
            if (i.shadow)
                for (u = 1; u <= i.lines; u++) p(u, -2, "progid:DXImageTransform.Microsoft.Blur(pixelradius=2,makeshadow=1,shadowopacity=.3)");
            for (u = 1; u <= i.lines; u++) p(u);
            return n(e, l)
        };
        c.prototype.opacity = function(t, e, i, o) {
            var n = t.firstChild;
            o = o.shadow && o.lines || 0;
            if (n && e + o < n.childNodes.length) {
                n = n.childNodes[e + o];
                n = n && n.firstChild;
                n = n && n.firstChild;
                if (n) n.opacity = i
            }
        }
    }
    var m = f(o("group"), { behavior: "url(#default#VML)" });
    if (!a(m, "transform") && m.adj) h();
    else i = a(m, "animation");
    return c
});
