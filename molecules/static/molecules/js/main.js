function getCookie(name) {
  let cookieValue = null;
  if (document.cookie && document.cookie !== '') {
    const cookies = document.cookie.split(';');
    for (let i = 0; i < cookies.length; i++) {
      const cookie = cookies[i].trim();
      if (cookie.substring(0, name.length + 1) === (name + '=')) {
        cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
        break;
      }
    }
  }
  return cookieValue;
}

function generateExpires(days) {
  const date = new Date();
  date.setTime(date.getTime() + days * 24 * 60 * 60 * 1000);
  return "; expires=" + date.toUTCString();
}

function setCookie(name, smiles, days) {
  const expires = days ? generateExpires(days) : "";

  if (name && smiles) {
    alert('Введіть або назву речовини, або smiles, але не обидва!');
    return;
  }

  if (!name && smiles) {
    CheckSmilesCode(smiles, expires);
  } else if (name && !smiles) {
    document.cookie = `moleculeName=${name}${expires}; path=/`;
    const link = document.getElementById("moleculeLinkName");
    link.href = window.location.origin + '/view_molecule/';
    link.style.display = "flex";
  } else {
    alert('Введіть хоча б одне значення — назву речовини або smiles!');
  }
}

function SendNameAndSmiles() {
  const nameInput = document.getElementById("inputName");
  const smilesInput = document.getElementById("inputSmiles");

  const name = nameInput.value.trim();
  const smiles = smilesInput.value.trim();

  if (!name && !smiles) {
    console.error("Помилка: поле введення порожнє!");
    alert("Будь ласка, введіть коректно назву молекули або smiles-код!");
    return;
  }

  setCookie(name, smiles, 2);

  nameInput.value = "";
  smilesInput.value = "";
}
function CheckSmilesCode(smiles, expires) {
  if (!smiles) {
    alert("Будь ласка, введіть smiles-код молекули!");
    return;
  }

  fetch(window.location.origin + "/check-smiles-code/", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      'X-CSRFToken': getCookie('csrftoken'),
    },
    body: JSON.stringify({ smiles: smiles }),
  })
    .then((response) => {
      if (!response.ok) { // Якщо статус не ОК, викидаємо помилку
        throw new Error(`HTTP Error! Status: ${response.status}`);
      }
      return response.json(); // Отримуємо JSON-дані
    })
    .then((data) => {
      if (data.status === "success") {
        document.cookie = `moleculeSmiles=${smiles}${expires}; path=/`;
        const link = document.getElementById("moleculeLinkSmileCode");
        link.href = window.location.origin + '/view_molecule/';
        link.style.display = "flex";
      }
    })
    .catch((error) => {
      // Тут catch вже коректно обробляє будь-які помилки, які виникають у fetch
      console.error("Помилка при запиті:", error.message); // Виводимо помилку в консоль
      alert("Будь ласка, введіть правильний SMILES-код або зверніться до інших ресурсів для його правильного написання!"); // Виводимо повідомлення користувачу
    });
}

function OpenModel(name, smiles){
	const expires = 7 ? generateExpires(7) : "";
	document.cookie = `moleculeName=${name}${expires}; path=/`;
	document.cookie = `moleculeSmiles=${smiles}${expires}; path=/`;

	window.location.href = window.location.origin + '/view_molecule/'
}