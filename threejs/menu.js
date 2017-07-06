function CreateMenu() {
    var MenuItems = new function() {
        this.EnableMouseTrackBall = true;
        this.EnableMouseActions = true;
        this.RestoreMouseActions = function() {
            restoreMouseOperations()
        };
        this.VolumeLabel = true;

        this.MaxNumberOfVolumes = 0;
        this.MaxLevel = 0;

        this.loadFromJsonFile = function() {
            fileInput.click();

            //load();

        };
        this.loadClones = function() {

            loadClones(jsClonesList);

        };
        this.setSelection = "";
        this.drawSelected = function() {
            var vol = MenuItems.setSelection;
            if (vol == "") {
                root.traverse(function(e) {
                    if (e instanceof THREE.Mesh) {
                        e.material.transparent = false;
                        e.material.opacity = 0.6;
                    }
                });

            } else {

                root.traverse(function(e) {
                    if (e instanceof THREE.Mesh) {
                        e.material.transparent = true;
                        e.material.opacity = 0.1;
                    }
                });
                var vId = root.getObjectByName(vol);

                if (vId instanceof THREE.Mesh) {
                    var mat = vId.material.clone();
                    vId.material = mat;
                    vId.material.transparent = false;
                    vId.material.opacity = 0.6;
                } else if (vId instanceof THREE.Group) {

                }
            }
        };

        this.hideVolumes = function() {
            for (var i = 0; i < scenes[0].children.length; i++) {
                if (scenes[0].children[i].name == "root") {
                    scenes[0].remove(scenes[0].children[i]);
                    break;
                }
            }
        };
        this.showVolumes = function() {
            for (var i = 0; i < scenes[0].children.length; i++) {
                var inScene = false;
                if (scenes[0].children[i].name == "root") {
                    inScene = true;
                    break;
                }
                if (!inScene) scenes[0].add(root);
            }
        };



        this.ShowTree = false;
        this.compactTree = false;

        this.traceParticles = false;
        this.restartParticles = false;

        this.particlesSizeScaleFactor = 20.;
        this.ReactedParticlesSizeScaleFactor = 20.;
        this.particleSpeedScaleFactor = 1.;


        this.gamma = true;
        this.electron = true;
        this.positron = true;
        this.proton = true;
        this.neutron = true;
        this.muon = true;
        this.deuterium = true;
        this.alpha = true;
        this.ion = true;
        this.energyMinCut = 0.;
        this.particleTransparency = true;
        this.particleOpacity = 0.3;


        this.localClipEnabled = false;
        this.clipingPlaneX = 0.
        this.clipingPlaneY = 0.
        this.clipingPlaneZ = 0.
            //===========================================================================
        this.ambientLightColor = "#99ff99";
        this.backgroundColor = "#505050";

        this.spotLightColor = "#ffffff";
        this.spotAngle = Math.PI / 5;
        this.spotPenumbra = 0.2;
        this.spotPositionX = 12000;
        this.spotPositionY = 12000;
        this.spotPositionZ = -12000;

        this.spotCastShadow = true;
        this.spotShadowCameraNear = 3;
        this.spotShadowCameraFar = 10;
        this.spotShadowMapSizeWidth = 1024;
        this.spotShadowMapSizeHeight = 1024;

        this.dirLightColor = "#ffffff";
        this.dirLightIntensity = 1;
        this.dirPositionX = 0;
        this.dirPositionY = 3000;
        this.dirPositionZ = 0;
        this.dirCastShadow = true;
        this.dirShadowCameraNear = 1;
        this.dirShadowCameraFar = 10;
        this.dirShadowCameraRight = 1;
        this.dirShadowCameraLeft = -1;
        this.dirShadowCameraTop = 1;
        this.dirShadowCameraBottom = -1;
        this.dirShadowMapSizeWidth = 1024;
        this.dirShadowMapSizeHeight = 1024;
    };

    //===========================================================================

    var gui = new dat.GUI({
        autoplace: false,
        width: 500
    });
    var fMain = gui.addFolder('Main Selection');
    fMain.add(MenuItems, 'EnableMouseTrackBall').onChange(function(newValue) {
        trackballControls.enabled = newValue;
    });
    fMain.add(MenuItems, 'EnableMouseActions')
    fMain.add(MenuItems, 'RestoreMouseActions')
    fMain.add(MenuItems, 'VolumeLabel');

    fMain.add(MenuItems, 'MaxNumberOfVolumes', 0, 3000000).listen();
    fMain.add(MenuItems, 'MaxLevel', 0, 20).listen();
    fMain.add(MenuItems, 'loadFromJsonFile');
    fMain.add(MenuItems, 'loadClones');

    fMain.add(MenuItems, 'hideVolumes');
    fMain.add(MenuItems, 'showVolumes');
    fMain.open();

    var fPsearchVolume = gui.addFolder('SearchVolume');
    fPsearchVolume.add(MenuItems, 'setSelection').onChange(function(newValue) {
        MenuItems.setSelection = newValue;

    });
    fPsearchVolume.add(MenuItems, 'drawSelected');

    var fPtreeStructure = gui.addFolder('TreeStructure');

    fPtreeStructure.add(MenuItems, 'ShowTree').onChange(function(newValue) {
        MenuItems.ShowTree = newValue;
        var r = scenes[0].getObjectByName("root");
        var p = scenes[0].getObjectByName("ParticlesInScene");
        var t = scenes[0].getObjectByName("ParticlesTraceInScene");
        var s = scenes[0].getObjectByName("treeSprites");
        var l = scenes[0].getObjectByName("treeLines");

        if (MenuItems.ShowTree) {

            //if (!treeCreated) {
                createTreeVolumes(maxLevel);
            //}
            if (r) scenes[0].remove(r);
            if (p) scenes[0].remove(p);
            if (t) scenes[0].remove(t);
            scenes[0].add(treeLines);
            scenes[0].add(treeSprites);
        } else {
            if (s) scenes[0].remove(s);
            if (s) scenes[0].remove(l);
            if (!r) scenes[0].add(root);
            if (!p) scenes[0].add(ParticlesInScene);
            if (!t) scenes[0].add(ParticlesTraceInScene);
        }

    }).listen();

    fPtreeStructure.add(MenuItems, 'compactTree').onChange(function(newValue) {
        MenuItems.compactTree = newValue;

        var r = scenes[0].getObjectByName("root");
        var p = scenes[0].getObjectByName("ParticlesInScene");
        var t = scenes[0].getObjectByName("ParticlesTraceInScene");
        var s = scenes[0].getObjectByName("treeSprites");
        var l = scenes[0].getObjectByName("treeLines");


        if (MenuItems.compactTree)    MenuItems.ShowTree = newValue;


        if (MenuItems.ShowTree ) {

            createTreeVolumes(maxLevel);
            if (r) scenes[0].remove(r);
            if (p) scenes[0].remove(p);
            if (t) scenes[0].remove(t);
            scenes[0].add(treeLines);
            scenes[0].add(treeSprites);
        } else {
            if (s) scenes[0].remove(s);
            if (s) scenes[0].remove(l);
            if (!r) scenes[0].add(root);
            if (!p) scenes[0].add(ParticlesInScene);
            if (!t) scenes[0].add(ParticlesTraceInScene);
        }
        


    });

    var fPselection = gui.addFolder('Particles');

    fPselection.add(MenuItems, 'traceParticles').onChange(function(newValue) {
        MenuItems.traceParticles = newValue;

    });
    fPselection.add(MenuItems, 'restartParticles').onChange(function(newValue) {
        MenuItems.restartParticles = newValue;
        if (newValue) RestartParticles();
    });
    fPselection.add(MenuItems, 'particlesSizeScaleFactor', {
        VerySmall: 0.1,
        Small: 0.2,
        HalfSize: 0.5,
        Normal: 1.,
        DoubleSize: 2.,
        Big: 5.,
        VeryBig: 10.,
        JUMBO: 20.,
        Giant: 50.
    });
    fPselection.add(MenuItems, 'ReactedParticlesSizeScaleFactor', {
        VerySmall: 0.1,
        Small: 0.2,
        HalfSize: 0.5,
        Normal: 1.,
        DoubleSize: 2.,
        Big: 5.,
        VeryBig: 10.,
        JUMBO: 20.,
        Giant: 50.
    });
    fPselection.add(MenuItems, 'particleSpeedScaleFactor', {
        Stopped: 0.,
        Boring: 0.01,
        microView: 0.05,
        closeView: 0.1,
        verySlow: 0.3,
        Slower: 0.5,
        Slow: 0.8,
        Normal: 1.
    });



    fPselection.add(MenuItems, 'gamma').onChange(function(newValue) {
        ParticleKindVisible[1] = newValue
    });
    fPselection.add(MenuItems, 'electron').onChange(function(newValue) {
        ParticleKindVisible[2] = newValue
    });
    fPselection.add(MenuItems, 'positron').onChange(function(newValue) {
        ParticleKindVisible[3] = newValue
    });
    fPselection.add(MenuItems, 'proton').onChange(function(newValue) {
        ParticleKindVisible[4] = newValue
    });
    fPselection.add(MenuItems, 'neutron').onChange(function(newValue) {
        ParticleKindVisible[5] = newValue
    });
    fPselection.add(MenuItems, 'muon').onChange(function(newValue) {
        ParticleKindVisible[6] = newValue
    });
    fPselection.add(MenuItems, 'deuterium').onChange(function(newValue) {
        ParticleKindVisible[7] = newValue
    });
    fPselection.add(MenuItems, 'alpha').onChange(function(newValue) {
        ParticleKindVisible[8] = newValue
    });
    fPselection.add(MenuItems, 'ion').onChange(function(newValue) {
        ParticleKindVisible[9] = newValue
    });


    var fTuning = gui.addFolder('Tuning Parameters');
    fTuning.add(MenuItems, 'energyMinCut', 0., 1.);
    fTuning.add(MenuItems, 'particleTransparency');
    fTuning.add(MenuItems, 'particleOpacity', 0., 1.);

    var fLocalClip = gui.addFolder("Local Clipping");

    fLocalClip.add(MenuItems, 'localClipEnabled').onChange(function(newValue) {
        MenuItems.localClipEnabled = newValue;
        renderer.localClippingEnabled = newValue;
    });
    fLocalClip.add(MenuItems, 'clipingPlaneX', -5000, 5000).onChange(function(newValue) {
        MenuItems.clipingPlaneX = newValue;
        clipPlanes[0].constant = newValue;
    });
    fLocalClip.add(MenuItems, 'clipingPlaneY', -5000, 5000).onChange(function(newValue) {
        MenuItems.clipingPlaneY = newValue;
        clipPlanes[1].constant = newValue;
    });
    fLocalClip.add(MenuItems, 'clipingPlaneZ', -5000, 5000).onChange(function(newValue) {
        MenuItems.clipingPlaneZ = newValue;
        clipPlanes[2].constant = newValue;
    });


    //===========================================================================
    var fAmbientLighting = gui.addFolder("Ambient Light");
    fAmbientLighting.addColor(MenuItems, 'ambientLightColor').onChange(function(newValue) {
        ambientLight.color.set(newValue);

    });
    fAmbientLighting.addColor(MenuItems, 'backgroundColor').onChange(function(newValue) {
        renderer.setClearColor(new THREE.Color(newValue));

    });

    var fSpotLighting = gui.addFolder("Spot Light");

    fSpotLighting.addColor(MenuItems, 'spotLightColor').onChange(function(newValue) {
        MenuItems.spotLightColor = newValue;
        spotLight.color.set(newValue);;
    });
    fSpotLighting.add(MenuItems, 'spotAngle').onChange(function(newValue) {
        MenuItems.spotAngle = newValue;
        spotLight.spotAngle = newValue;
    });
    fSpotLighting.add(MenuItems, 'spotPenumbra').onChange(function(newValue) {
        MenuItems.spotPenumbra = newValue;
        spotLight.spotPenumbra = newValue;
    });
    fSpotLighting.add(MenuItems, 'spotPositionX').onChange(function(newValue) {
        MenuItems.spotPositionX = newValue;
        spotLight.position.x = newValue;
    });
    fSpotLighting.add(MenuItems, 'spotPositionY').onChange(function(newValue) {
        MenuItems.spotPositionY = newValue;
        spotLight.position.y = newValue;
    });
    fSpotLighting.add(MenuItems, 'spotPositionZ').onChange(function(newValue) {
        MenuItems.spotPositionZ = newValue;
        spotLight.position.z = newValue;
    });

    fSpotLighting.add(MenuItems, 'spotCastShadow').onChange(function(newValue) {
        MenuItems.spotCastShadow = newValue;
        spotLight.castShadow = newValue;
    });
    fSpotLighting.add(MenuItems, 'spotShadowCameraNear').onChange(function(newValue) {
        MenuItems.spotShadowCameraNear = newValue;
        spotLight.shadow.camera.near = newValue;
    });
    fSpotLighting.add(MenuItems, 'spotShadowCameraFar').onChange(function(newValue) {
        MenuItems.spotShadowCameraFar = newValue;
        spotLight.shadow.camera.far = newValue;
    });
    fSpotLighting.add(MenuItems, 'spotShadowMapSizeWidth').onChange(function(newValue) {
        MenuItems.spotShadowMapSizeWidth = newValue;
        spotLight.shadow.mapSize.width = newValue;
    });
    fSpotLighting.add(MenuItems, 'spotShadowMapSizeHeight').onChange(function(newValue) {
        MenuItems.spotShadowMapSizeWidth = newValue;
        spotLight.shadow.mapSize.height = newValue;
    });

    var fDirLighting = gui.addFolder("Diectional Light");

    fDirLighting.addColor(MenuItems, 'dirLightColor').onChange(function(newValue) {
        MenuItems.dirLightColor = newValue;
        dirLight.color.set(newValue);
    });
    fDirLighting.add(MenuItems, 'dirLightIntensity').onChange(function(newValue) {
        MenuItems.dirLightIntensity = newValue;
        dirLight.intensity = newValue;
    });
    fDirLighting.add(MenuItems, 'dirPositionX').onChange(function(newValue) {
        MenuItems.dirPositionX = newValue;
        dirLight.position.x = newValue;
    });
    fDirLighting.add(MenuItems, 'dirPositionY').onChange(function(newValue) {
        MenuItems.dirPositionY = newValue;
        dirLight.position.y = newValue;
    });
    fDirLighting.add(MenuItems, 'dirPositionZ').onChange(function(newValue) {
        MenuItems.dirPositionZ = newValue;
        dirLight.position.z = newValue;
    });
    fDirLighting.add(MenuItems, 'dirCastShadow').onChange(function(newValue) {
        MenuItems.dirCastShadow = newValue;
        dirLight.castShadow = newValue;
    });
    fDirLighting.add(MenuItems, 'dirShadowCameraNear').onChange(function(newValue) {
        MenuItems.dirCastShadow = newValue;
        dirLight.castShadow = newValue;
    });
    fDirLighting.add(MenuItems, 'dirShadowCameraFar').onChange(function(newValue) {
        MenuItems.dirShadowCameraFar = newValue;
        dirLight.shadow.camera.far = newValue;
    });
    fDirLighting.add(MenuItems, 'dirShadowCameraRight').onChange(function(newValue) {
        MenuItems.dirShadowCameraRight = newValue;
        dirLight.shadow.camera.right = newValue;
    });
    fDirLighting.add(MenuItems, 'dirShadowCameraLeft').onChange(function(newValue) {
        MenuItems.dirShadowCameraLeft = newValue;
        dirLight.shadow.camera.left = newValue;
    });
    fDirLighting.add(MenuItems, 'dirShadowCameraTop').onChange(function(newValue) {
        MenuItems.dirShadowCameraTop = newValue;
        dirLight.shadow.camera.top = newValue;
    });
    fDirLighting.add(MenuItems, 'dirShadowCameraBottom').onChange(function(newValue) {
        MenuItems.dirShadowCameraBottom = newValue;
        dirLight.shadow.camera.bottom = newValue;
    });
    fDirLighting.add(MenuItems, 'dirShadowMapSizeWidth').onChange(function(newValue) {
        MenuItems.dirShadowMapSizeWidth = newValue;
        dirLight.shadow.mapSize.width = newValue;
    });
    fDirLighting.add(MenuItems, 'dirShadowMapSizeHeight').onChange(function(newValue) {
        MenuItems.dirShadowMapSizeHeight = newValue;
        dirLight.shadow.mapSize.height = newValue;
    });

    //===========================================================================



    //f2.open();
    //gui.closed = true;
    return (MenuItems);
}
