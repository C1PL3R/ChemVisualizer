const molecule = document.getElementById("molecule").innerText;
const molecule_name = document.getElementById("molecule_name").innerText;

document.addEventListener("DOMContentLoaded", function () {
    var molBlock = molecule;
    if (molBlock.trim()) {
        var viewer = $3Dmol.createViewer("viewer", { backgroundColor: "transparent" });
        viewer.addModel(molBlock, "sdf");

        // Тонкі зв'язки та невеликі сфери
        viewer.setStyle({}, { stick: { radius: 0.07 }, sphere: { radius: 0.3 } });

        viewer.setBackgroundColor("white");


        viewer.zoomTo();
        viewer.render();
        viewer.zoom(1.2, 500);

        document.getElementById('animation').addEventListener('change', function () {
            var text = document.getElementById('turn-on-off');
            if (this.checked) {
                viewer.spin(true, 2);
                text.innerHTML = "Вимкнути анімацію";
            } else {
                viewer.spin(false);
                text.innerHTML = "Увімкнути анімацію";
            }
        });

        document.getElementById('show-text').addEventListener('change', function () {
            var text = document.getElementById('show-text-span');

            if (this.checked) {
                var atoms = viewer.getModel().selectedAtoms({});
                atoms.forEach(atom => {
                    viewer.addLabel(atom.elem, {
                        position: { x: atom.x, y: atom.y, z: atom.z },
                        fontSize: 18,
                        fontColor: "black",
                        backgroundColor: "white",  // Для перевірки можна поставити колір
                        backgroundOpacity: 0.0      // Повна прозорість фону
                    });
                });

                text.innerHTML = "Приховати підписи";
            } else {
                viewer.removeAllLabels();
                text.innerHTML = "Показати підписи"; 
            }
            viewer.render();
        });

    } else {
        console.error("MolBlock порожній або не був переданий.");
    }
});

document.getElementById("downloadBtn").addEventListener("click", function () {
    var Data = molecule;
    if (Data.trim()) {
        var blob = new Blob([Data], { type: "text/plain" });
        var link = document.createElement("a");
        link.href = URL.createObjectURL(blob);
        link.download = molecule_name + ".sdf";
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    } else {
        console.error("SDF дані порожні!");
    }
});


function openFullscreen() {
    const elem = document.getElementById("viewer");
    // Відкрити fullscreen
    if (elem.requestFullscreen) {
        elem.requestFullscreen();
    } else if (elem.webkitRequestFullscreen) {
        elem.webkitRequestFullscreen();
    } else if (elem.msRequestFullscreen) {
        elem.msRequestFullscreen();
    }
}

