function dc = chargerRF(discount2,CC)

NECL = 10*365;

id = discount2 / (365);

dc = (id*(CC*(1+id)^NECL))/ ( (1+id)^NECL - 1);

end