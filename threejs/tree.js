function leafPosSet(node, position) {
    //if (node.compacted === undefined) 
    node.compacted = false;
    for (var i = 0; i < node.children.length; i++) {
        position = leafPosSet(node.children[i], position);
    }
    var pos = 1;
    var cmpt= 0;
    if (node.children.length == 0) {

        if (MenuItems.compactTree) {

            if (!node.parent.compacted) {

                if (node.parent.children.length > treeNodesToCompact) {
                    var n = 0;
                    for (var i = 0; i < node.parent.children.length; i++) {
                        if (node.parent.children[i].children.length == 0) n++;
                        if (n > treeNodesToCompact) {
                            node.parent.compacted = true;
                            pos = 1;
                            break;
                        }
                    }
                }
            } else {
                node.compacted = true;
                pos = 0;
                cmpt= -1;
            }

        }
        node.treePosition = position+cmpt;

        position += pos;

        //console.log("Node: ", node.name, "   no of Childs: ", node.children.length, " Tree position", node.treePosition, " next pos:", position);
    }
    return position;
}

function nodesPosSet(node) {
    for (var i = 0; i < node.children.length; i++) {
        nodesPosSet(node.children[i]);
    }
    if (node.children.length > 0) {

        if (node.children.length == 1) {
            node.treePosition = node.children[0].treePosition;
        } else {
            node.treePosition = (node.children[0].treePosition + node.children[node.children.length - 1].treePosition) / 2.;
        }
    }

}

function CreateNodesSprites(node, startNode, maxPosition, maxLevel, level) {

    var label = node.name;
    var color=treeSpriteColors[0];
    if (MenuItems.compactTree) {
        if (node != startNode) {
            if (node.parent.compacted && node.children.length == 0) {
                label = "Compact...";
                var color=treeSpriteColors[1];
            }
        }
    }
//    var newSprite = makeTextSprite(label, { fontsize: 24, backgroundColor: { r: 0, g: 100, b: 255, a: 1 } });

    var newSprite = makeTextSprite(label, color);
    newSprite.nodeId = node;

    newSprite.position.x = 12 * (-1 * maxPosition / 2 + node.treePosition);
    newSprite.position.y = 20 * (-1 * maxLevel / 2 + level);
    newSprite.position.z = 0.;
    //console.log(newSprite.position.x, newSprite.position.y, newSprite.position.z);
    if (node.children.length > 0) {

        var vLinePosX1 = 12 * (-1 * maxPosition / 2 + node.children[0].treePosition);
        var vLinePosX2 = 12 * (-1 * maxPosition / 2 + node.children[node.children.length - 1].treePosition);
        var vLinePosY = 20 * (-1 * maxLevel / 2 + level + 0.5)

        //vertical line for the childrens of this node (up)
        TreeLineGeometry.vertices.push(new THREE.Vector3(newSprite.position.x, newSprite.position.y, 0.));
        TreeLineGeometry.vertices.push(new THREE.Vector3(newSprite.position.x, vLinePosY, 0.));

        if (node.children.length > 1) {
            //horizontal line for the childrens of this node (up)
            TreeLineGeometry.vertices.push(new THREE.Vector3(vLinePosX1, vLinePosY, 0.));
            TreeLineGeometry.vertices.push(new THREE.Vector3(vLinePosX2, vLinePosY, 0.));
            //console.log("  Line up to", vLinePosY, " + Horizontal line: From:", vLinePosX1, vLinePosY, 0., " To: ", vLinePosX2, vLinePosY, 0.);

        } else {
            //console.log("  Line up to", vLinePosY);

        }

    }

    if (node != startNode) {
        var vLinePosYn = 20 * (-1 * maxLevel / 2 + level - 0.5);

        TreeLineGeometry.vertices.push(new THREE.Vector3(newSprite.position.x, newSprite.position.y, 0.));

        TreeLineGeometry.vertices.push(new THREE.Vector3(newSprite.position.x, vLinePosYn, 0.));
        //console.log("  Line down to:", vLinePosYn);
    }
    //console.log("-----------------Next node--------------");

    treeSprites.add(newSprite);

    for (var i = 0; i < node.children.length; i++) {
        if (MenuItems.compactTree) {
            if (node.children[i].compacted && node.children[i].children.length == 0) {
                //console.log("Node compacted: ", node.name)
            } else {
                CreateNodesSprites(node.children[i], startNode, maxPosition, maxLevel, level + 1);
            }
        } else {
            CreateNodesSprites(node.children[i], startNode, maxPosition, maxLevel, level + 1);
        }
    }
}




function makeTextSprite(message, bgcolor, parameters) {
    if (parameters === undefined) parameters = {};

    var fontface = parameters.hasOwnProperty("fontface") ?
        parameters["fontface"] : "Arial";

    var fontsize = parameters.hasOwnProperty("fontsize") ?
        parameters["fontsize"] : 18;

    var borderThickness = parameters.hasOwnProperty("borderThickness") ?
        parameters["borderThickness"] : 4;

    var borderColor = parameters.hasOwnProperty("borderColor") ?
        parameters["borderColor"] : { r: 0, g: 0, b: 0, a: 1.0 };

    var backgroundColor = parameters.hasOwnProperty("backgroundColor") ?
        parameters["backgroundColor"] : { r: 255*bgcolor.r, g: 255*bgcolor.g, b: 255*bgcolor.b, a: 1.0 };

    var canvasTree = document.createElement('canvas');
    canvasTree.width = 128.;
    canvasTree.height = 64.;

    var context = canvasTree.getContext('2d');

    context.font = "Bold " + fontsize + "px " + fontface;

    // get size data (height depends only on font size)
    var metrics = context.measureText(message);
    var textWidth = metrics.width;

    // background color
    context.fillStyle = "rgba(" + backgroundColor.r + "," + backgroundColor.g + "," + backgroundColor.b + "," + backgroundColor.a + ")";
    // border color
    context.strokeStyle = "rgba(" + borderColor.r + "," + borderColor.g + "," + borderColor.b + "," + borderColor.a + ")";

    context.lineWidth = borderThickness;
    //roundRect(context, borderThickness/2, borderThickness/2, textWidth + borderThickness, fontsize * 1.4 + borderThickness, 6);
    roundRect(context, borderThickness / 2, borderThickness / 2, 128 - borderThickness, 64 - borderThickness, 6);
    // 1.4 is extra height factor for text below baseline: g,j,p,q.

    // text color
    context.fillStyle = "rgba(0, 0, 0, 1.0)";

    context.fillText(message, borderThickness, fontsize + borderThickness);

    // canvasTree contents will be used for a texture
    var textureTree = new THREE.Texture(canvasTree);
    textureTree.needsUpdate = true;

    var spriteMaterial = new THREE.SpriteMaterial({ map: textureTree });
    var sprite = new THREE.Sprite(spriteMaterial);
    sprite.scale.set(10, 5, 1.0);
    return sprite;
}

// function for drawing rounded rectangles
function roundRect(ctx, x, y, w, h, r) {
    //console.log( "CANVAS: ",x, y, w, h, r);
    ctx.beginPath();
    ctx.moveTo(x + r, y);
    ctx.lineTo(x + w - r, y);
    ctx.quadraticCurveTo(x + w, y, x + w, y + r);
    ctx.lineTo(x + w, y + h - r);
    ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
    ctx.lineTo(x + r, y + h);
    ctx.quadraticCurveTo(x, y + h, x, y + h - r);
    ctx.lineTo(x, y + r);
    ctx.quadraticCurveTo(x, y, x + r, y);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
}

function disposeGeoTree() {
    scenes[0].remove(treeLines);
    scenes[0].remove(treeSprites);
    for (var i = treeLines.children.length - 1; i >= 0; i--) {
        var g = treeLines.children[i].geometry;
        var m = treeLines.children[i].material;

        treeLines.remove(treeLines.children[i]);
        g.dispose();
        m.dispose();

    }
    TreeLineGeometry.dispose();
    TreeLineGeometry = new THREE.Geometry();


    for (var i = treeSprites.children.length - 1; i >= 0; i--) {
        //var g = treeSprites.children[i].geometry;
        var m = treeSprites.children[i].material;
        var t = treeSprites.children[i].material.map; //.image;
        treeSprites.children[i].material.map = null;
        treeSprites.remove(treeSprites.children[i]);
        //g.dispose();
        t.dispose();
        m.dispose();
    }
}
//geometry.dispose();
//material.dispose();
//texture.dispose();

function createTreeVolumes(maxLevel) {
    if (treeCreated) {
        disposeGeoTree();
    }
    treeCompact = MenuItems.compactTree;
    treeCreated = true;
    var position = 0;
    var maxlev = 0;
    var level = 0;
    //count the number of leafs = node.treePosition
    console.log("Leafs position .... ");
    position = leafPosSet(root, position);
    console.log("Total number of leaf nodes: ", position);
    //calculate the position of mother nodes (middle distance between 1st and last son)
    console.log("prosition nodes .... ");
    nodesPosSet(root);
    console.log("prosition nodes DONE: ");

    CreateNodesSprites(root, root, position, maxLevel, level);
    TreeLineGeometry.computeLineDistances();

    var TreeLineMaterial = new THREE.LineDashedMaterial({ color: 0xffaa00, dashSize: 1., gapSize: 0.5 });
    var TreeLine = new THREE.LineSegments(TreeLineGeometry, TreeLineMaterial);
    treeLines.add(TreeLine);

}

/*function disposeGeoTree() {
    for (var i = treeLines.children.length - 1; i >= 0; i--) {
        var g = treeLines.children[i].geometry;
        var m = treeLines.children[i].material;

        scenes[0].remove(treeLines.children[i]);
        treeLines.remove(treeLines.children[i]);
        g.dispose();
        m.dispose();
    }

    for (var i = treeSprites.children.length - 1; i >= 0; i--) {
        //var g = treeSprites.children[i].geometry;
        var m = treeSprites.children[i].material;
        var t = treeSprites.children[i].material.map;
        treeSprites.children[i].material.map = null;
        scenes[0].remove(treeSprites.children[i]);
        treeSprites.remove(treeSprites.children[i]);
        //g.dispose();
        m.dispose();
        t.dispose();
    }
}
*/
