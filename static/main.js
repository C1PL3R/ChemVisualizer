function sendMoleculeRequest(smiles) {
  fetch("http://192.168.0.108:8000/molecule-smiles/", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({ name: smiles.trim() }), // Виправлено ключ
  })
    .then((response) => {
      if (!response.ok) {
        throw new Error(`HTTP Error! Status: ${response.status}`);
      }
      return response.json();
    })
    .then((data) => {
      if (data.status === "success") {
        console.log(data.smiles);
        var link = document.getElementById("moleculeLink");

        link.href = `view_molecule/${data.name}/${data.smiles}`;
        link.style.display = "flex";

      } else {
        console.error(`Помилка сервера: ${data.error}`);
      }
    })
    .catch((error) => {
      console.error("Помилка запиту:", error.message);
      alert("Введіть назву речовини коректно!"); // Виведе сповіщення
    });
}

function SendName() {
  let name = document.getElementById("inputName").value.trim();

  if (!name) {
    console.error("Помилка: поле введення порожнє!");
    alert("Будь ласка, введіть назву молекули!");
    return;
  }

  sendMoleculeRequest(encodeURIComponent(name)); // Кодуємо спеціальні символи
}
