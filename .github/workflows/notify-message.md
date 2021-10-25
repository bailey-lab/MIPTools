Please build the container and deploy it to the Sylabs Cloud. It has been two weeks since the last reminder.

Steps to build:

```bash
sudo singularity build miptools.sif MIPTools.def
```

Steps to deploy:

```bash
singularity remote login
singularity sign miptools.sif
singularity push miptools.sif library://apascha1/MIPTools/miptools:{tag}
```
