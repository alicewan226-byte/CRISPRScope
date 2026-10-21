This [Code Ocean](https://codeocean.com) Compute Capsule will allow you to run and reproduce the results of [CRISPRScope – single-cell CRISPR perturbation analysis](https://codeocean.ir-scc-codeocean-prod.awscloud.abbvienet.com/capsule/6234214/tree) on your local machine<sup>1</sup>. Follow the instructions below, or consult [our knowledge base](https://docs.codeocean.com/user-guide/compute-capsule-basics/managing-capsules/exporting-capsules-to-your-local-machine) for more information. Don't hesitate to reach out to [Support](mailto:support@codeocean.com) if you have any questions.

<sup>1</sup> You may need access to additional hardware and/or software licenses.

# Prerequisites

- [Docker Community Edition (CE)](https://www.docker.com/community-edition)

# Instructions

## Log in to the Docker registry

In your terminal, execute the following command, providing your password or API key when prompted for it:
```shell
docker login -u alice.wan@abbvie.com registry.codeocean.ir-scc-codeocean-prod.awscloud.abbvienet.com
```

## Run the Capsule to reproduce the results

In your terminal, navigate to the folder where you've extracted the Capsule and execute the following command, adjusting parameters as needed:
```shell
docker run --platform linux/amd64 --rm \
  --workdir /code \
  --volume "$PWD/code":/code \
  --volume "$PWD/data":/data \
  --volume "$PWD/results":/results \
  registry.codeocean.ir-scc-codeocean-prod.awscloud.abbvienet.com/capsule/15f7d16e-6969-4933-a80f-ef84227db048 \
  bash run
```
