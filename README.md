# Stéganographie

## Énoncé

> Ça y est, Adèle, je viens de déposer mon brevet de Cryptostéganocodes.
> > Jamais entendu parler, ça sert à quoi ?

> C'est de la stéganographie, je mets l'information secrète au dessus d'un code
> correcteur, sous forme d'erreurs sciemment ajoutées.
> Personne ne soupçonne son existence.
> > Et ça fonctionne bien ?

> Vois plutôt, l'entier suivant représente un mot de code BCH 35-correcteur, de longueur 127 sur F128,
> je l'ai utilisé pour cacher des informations précieuses.
> 0xdeb4897ea42c90ed3d62dda7331153957e97ef93f8341f40b3cc37cbbd6b6a4f25ee7c1639723e2e9d8670c96599e48fc0c6c1cc5217f9d2f5c9db8442dfb5074fccc37ce1ba2f8a549b05d78fccb6de220216b62b776b4868e38c4864be26a5416ab0588b171bb613d1296d6b18
> > C'est quoi comme info ?

> Juste une chaîne ascii.

Quelques minutes plus tard...

> > Espèce d'idiot ! laisse tomber, tu n'as aucune chance.

## Format

Vous devez écrire un programme `main.gp` en Pari/GP dont l'exécution via
```
gp -fq < main.gp
```
affiche (uniquement) la solution cherchée.

taper `make check` permet de vérifier que votre solution est bonne.

Veillez à mettre des commentaires sur votre démarche et sa validité
dans le fichier ``main.gp``.

De manière alternative, vous pouvez écrire un programme `main.c` qui
fasse la même chose.

## Validation

Il reste à faire un commit et à envoyer votre solution pour l'enregistrer
```
git add main.gp
git commit -m 'ma solution'
git push
```

