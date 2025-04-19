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

function setMoleculeCookies(name, smiles) {
    const expiry = new Date();
    expiry.setTime(expiry.getTime() + (7 * 24 * 60 * 60 * 1000));
    const expires = `; path=/; expires=${expiry.toUTCString()}`;
    document.cookie = `moleculeSmiles=${smiles}${expires}`;
    document.cookie = `moleculeName=${name}${expires}`;
}

function OpenModel(name, smiles){
	setMoleculeCookies(name, smiles)

	window.location.href = window.location.origin + '/view_molecule/'
}

document.getElementById("sendNameBtn").addEventListener("click", function () {
	const inputName = document.getElementById("inputName");
	const loader = document.getElementById("nameLoader");
	const loader_container = document.getElementById("loader_container");

	var name = inputName.value.trim();
	if (!name) return;

	loader.style.display = "inline-block";
	loader_container.style.display = "flex"; // Показати завантаження

	fetch(window.location.origin + "/send-name/", {
		method: "POST",
		headers: {
			"Content-Type": "application/json",
			"X-CSRFToken": getCookie("csrftoken"),
		},
		body: JSON.stringify({ 'name': name }),
	})
	.then(response => response.json())
	.then((data) => {
		if (data.status === "success") {
			setMoleculeCookies( data.name, data.smiles);
			const sendBtn = document.getElementById("sendNameBtn");
			sendBtn.style.display = 'none';

			const link = document.getElementById("moleculeLinkName");
			link.href = window.location.origin + "/view_molecule/";
			link.style.display = "inline";
		}
	})
	.catch((error) => {
		console.error("Помилка при запиті:", error.message);
		alert("Будь ласка, введіть правильну назву або зверніться до інших ресурсів для її правильного написання!");
	})
	.finally(() => {
		loader.style.display = "none";
		loader_container.style.display = "none";
	});
});


document.getElementById("sendFormulaBtn").addEventListener("click", function () {
	const inputFormula = document.getElementById("inputFormula");
	const loader = document.getElementById("formulaLoader");
	const loader_container = document.getElementById("loader_container");

	var formula = inputFormula.value.trim();
	if (!formula) return;

	loader.style.display = "inline-block";
	loader_container.style.display = "flex"; // Показати завантаження

	fetch(window.location.origin + "/send-formula/", {
		method: "POST",
		headers: {
			"Content-Type": "application/json",
			"X-CSRFToken": getCookie("csrftoken"),
		},
		body: JSON.stringify({ 'formula': formula }),
	})
	.then(response => response.json())
	.then((data) => {
		if (data.status === "success") {
			setMoleculeCookies(data.name, data.smiles);
			const sendBtn = document.getElementById("sendFormulaBtn");
			sendBtn.style.display = 'none';

			const link = document.getElementById("moleculeLinkFormula");
			link.href = window.location.origin + "/view_molecule/";
			link.style.display = "inline";
		}
	})
	.catch((error) => {
		console.error("Помилка при запиті:", error);
		alert("Будь ласка, введіть правильну формулу або зверніться до інших ресурсів для її правильного написання!");
	})
	.finally(() => {
		loader.style.display = "none";
		loader_container.style.display = "none";
	});
});

var value = 'block';

document.getElementById("nameOfSubstanceBtn").addEventListener("click", function () {
	const nameOfSubstance = document.getElementById("nameOfSubstance");
	const formulaOfSubstance = document.getElementById("formulaOfSubstance");
	const input_data_container = document.getElementById("input_data_container");


	if (nameOfSubstance.style.display === value) {
		nameOfSubstance.style.display = "none";
		input_data_container.style.display = "none";
	} else {
		nameOfSubstance.style.display = value;
		input_data_container.style.display = value;
		formulaOfSubstance.style.display = "none";
	}
});

document.getElementById("formulaOfSubstanceBtn").addEventListener("click", function () {
	const formulaOfSubstance = document.getElementById("formulaOfSubstance");
	const nameOfSubstance = document.getElementById("nameOfSubstance");
	const input_data_container = document.getElementById("input_data_container");

	if (formulaOfSubstance.style.display === value) {
		formulaOfSubstance.style.display = "none";
		input_data_container.style.display = "none";
	} else {
		formulaOfSubstance.style.display = value;
		input_data_container.style.display = value;
		nameOfSubstance.style.display = "none";
	}
});

